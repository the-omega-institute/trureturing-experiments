#!/usr/bin/env python3
"""Bounded, synthetic, same-filesystem SC1 removal and interruption verification.

Run with --source pointing to the saved implementation worktree. Only native
ls-tree metadata is read from the pinned manifest. Fixture repositories contain
empty dummy commits, never a checkout/copy of project source or cache data.
"""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import stat
import struct
import subprocess
import sys
import tempfile
import time
from unittest.mock import patch

MANIFEST = "1533c8bdaf8c36bb1d6a28129a52e8e14d29c7b7"
LIMIT = 2 * 1024**3
BLOCK = 4096
SOURCE_PATHS = ["tools/scripts/worktree/" + n + ".py" for n in
                ("worktree_protocol", "worktree_pool", "worktree_preservation")]


def git(root, *args, env=None):
    result = subprocess.run(["git", "-C", str(root), *map(str, args)],
                            env={**os.environ, **(env or {})}, capture_output=True, check=True)
    return result.stdout.decode().strip()


def bindings(source):
    return dict(head=git(source, "rev-parse", "HEAD"),
                files={p: hashlib.sha256((source / p).read_bytes()).hexdigest() for p in SOURCE_PATHS})


def refs(main):
    return git(main, "for-each-ref", "--format=%(refname) %(objectname)", "refs/heads/")


def census(root):
    files = directories = logical = allocated = links = 0
    for parent, dirs, names in os.walk(root, followlinks=False):
        directories += 1
        for name in names + dirs:
            p = Path(parent) / name
            st = p.lstat()
            if stat.S_ISREG(st.st_mode):
                files += 1
                logical += st.st_size
                allocated += st.st_blocks * 512
            elif stat.S_ISLNK(st.st_mode):
                links += 1
    return dict(regular_files=files, directories=directories, logical_bytes=logical,
                allocated_regular_bytes=allocated, symlinks=links)


def state(protocol, main, tree):
    return dict(root=tree.exists(), pointer=(tree / ".git").exists(),
                registered=any(Path(row["worktree"]) == tree for row in protocol.inventory(main)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--parent", type=Path, default=Path("/Users/auricstudio"))
    args = parser.parse_args()
    source = args.source.resolve()
    assert sys.version_info >= (3, 12)
    assert platform.system() == "Darwin", "hole-punch and immutable-flag probe is macOS-specific"
    sys.path.insert(0, str(source / "tools/scripts/worktree"))
    import worktree_preservation as owner
    import worktree_protocol as protocol
    initial = bindings(source)
    raw = subprocess.check_output(["git", "-C", str(source), "ls-tree", "-rlz", MANIFEST])
    entries = []
    directories = set()
    fanout = {}
    for row in raw.split(b"\0"):
        if not row:
            continue
        fields, rawpath = row.split(b"\t", 1)
        mode, kind, oid, size = fields.split()
        path = Path(os.fsdecode(rawpath))
        assert not path.is_absolute() and ".." not in path.parts and path.parts[0] != ".git"
        entries.append((path, int(size), mode.decode()))
        directories.update(p for p in path.parents if str(p) != ".")
        fanout[str(path.parent)] = fanout.get(str(path.parent), 0) + 1
    assert len(entries) == 393567 and len(directories) == 6793
    assert max(fanout.values()) == 147197
    fixture = Path(tempfile.mkdtemp(prefix=".assignment5-sc1-", dir=args.parent)).resolve()
    assert fixture.stat().st_dev == source.stat().st_dev
    result = dict(source=initial, manifest=dict(revision=MANIFEST, sha256=hashlib.sha256(raw).hexdigest(),
                  files=len(entries), directories=len(directories), maximum_fanout=max(fanout.values()),
                  logical_bytes=sum(e[1] for e in entries), modes={m:sum(e[2] == m for e in entries)
                  for m in sorted(set(e[2] for e in entries))}),
                  environment=dict(python=sys.version, executable=sys.executable,
                  platform=platform.platform(), git=git(source, "--version"), device=fixture.stat().st_dev),
                  fixture=str(fixture), regular_allocation_limit=LIMIT, passed=False)
    original_run = subprocess.run
    calls = []

    def measured(command, *pos, **kw):
        watched = command[0] == "/bin/rm" or ("worktree" in command and "remove" in command)
        if not watched:
            return original_run(command, *pos, **kw)
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.monotonic()
        record = dict(phase="payload" if command[0] == "/bin/rm" else "final-native",
                      timeout_seconds=kw.get("timeout"))
        try:
            completed = original_run(command, *pos, **kw)
            record["exit"] = completed.returncode
            return completed
        except subprocess.CalledProcessError as error:
            record.update(exit=error.returncode, stderr=error.stderr.decode(errors="replace")[-4096:])
            raise
        finally:
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            record.update(wall_seconds=time.monotonic()-start,
                          user_cpu_seconds=after.ru_utime-before.ru_utime,
                          system_cpu_seconds=after.ru_stime-before.ru_stime)
            calls.append(record)

    try:
        dummy = fixture / "main"
        git(fixture, "init", "--initial-branch=dev", dummy)
        git(dummy, "config", "user.name", "SC1 disposable verification")
        git(dummy, "config", "user.email", "sc1@example.invalid")
        git(dummy, "commit", "--allow-empty", "-m", "dummy base")
        base = git(dummy, "rev-parse", "HEAD")
        old_env = dict(GIT_AUTHOR_DATE="@"+str(int(time.time())-7*3600),
                       GIT_COMMITTER_DATE="@"+str(int(time.time())-7*3600))

        def tree(name):
            target = fixture / name
            git(dummy, "worktree", "add", "-b", name, target, base)
            git(target, "commit", "--allow-empty", "-m", name, env=old_env)
            return target

        # Actual native final failure: immutable root pointer. No restoration or
        # registered-path fallback is attempted when native Git drops registration.
        final = tree("final-fault")
        (final / "payload").write_bytes(b"owned")
        before_refs = refs(dummy)
        pointer = final / ".git"
        os.chflags(pointer, stat.UF_IMMUTABLE)
        start_calls = len(calls)
        try:
            with patch.object(owner.subprocess, "run", side_effect=measured):
                try:
                    owner.remove_registered(dummy, final)
                    raise AssertionError("immutable pointer incorrectly succeeded")
                except protocol.GitFailure as error:
                    result["final_failure"] = dict(reason=str(error), state=state(protocol,dummy,final))
            assert result["final_failure"]["state"] == dict(root=True,pointer=True,registered=False)
            assert not (final / "payload").exists()
            try:
                owner.remove_registered(dummy, final)
                raise AssertionError("orphan was accepted")
            except protocol.Refused:
                pass
            assert refs(dummy) == before_refs
            result["final_failure"].update(native_calls=calls[start_calls:], refs_preserved=True,
                                           orphan_retry_refused=True)
        finally:
            os.chflags(pointer, 0)
            shutil.rmtree(final)  # This experiment owns every byte of this dummy orphan.

        interrupted = tree("interrupt-fault")
        payload = interrupted / "payload"
        payload.mkdir()
        for index in range(20000):
            (payload / str(index)).write_bytes(b"x")
        before_refs = refs(dummy)
        interrupt_record = {}

        def interrupt_rm(command, *pos, **kw):
            if command[0] != "/bin/rm":
                return original_run(command, *pos, **kw)
            before = resource.getrusage(resource.RUSAGE_CHILDREN)
            start = time.monotonic()
            # Exactly this owned child is signalled and synchronously reaped.
            with subprocess.Popen(command, cwd=kw["cwd"], env=kw["env"],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
                deadline = time.monotonic() + 10
                watched = [payload / str(index) for index in range(0,20000,79)]
                while all(p.exists() for p in watched):
                    assert time.monotonic() < deadline, "no confirmed payload progress"
                    try:
                        process.communicate(timeout=0.01)
                        raise AssertionError("interrupt fixture completed before signal")
                    except subprocess.TimeoutExpired:
                        pass
                process.send_signal(signal.SIGINT)
                out, err = process.communicate(timeout=10)
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                interrupt_record.update(exit=process.returncode, wall_seconds=time.monotonic()-start,
                    user_cpu_seconds=after.ru_utime-before.ru_utime,
                    system_cpu_seconds=after.ru_stime-before.ru_stime)
                assert process.returncode == -signal.SIGINT
                raise subprocess.CalledProcessError(process.returncode,command,output=out,stderr=err)

        with patch.object(owner.subprocess,"run",side_effect=interrupt_rm):
            try:
                owner.remove_registered(dummy,interrupted)
                raise AssertionError("interruption incorrectly succeeded")
            except protocol.GitFailure as error:
                interrupt_record["reported_reason"] = str(error)
        remaining = sum(1 for p in payload.iterdir())
        assert 0 < remaining < 20000
        protocol.identity(dummy,interrupted)
        interrupt_record.update(initial_files=20000,remaining_files=remaining,
                                state=state(protocol,dummy,interrupted))
        with patch.object(owner.subprocess,"run",side_effect=measured):
            owner.remove_registered(dummy,interrupted)
        assert not interrupted.exists() and refs(dummy)==before_refs
        interrupt_record.update(ordinary_qualified_retry=True,refs_preserved=True)
        result["interrupt"] = interrupt_record

        target = tree("representative")
        before_refs = refs(dummy)
        external = fixture / "external-sentinel"
        external.write_bytes(b"owned external sentinel\n")
        (target / ".lake").symlink_to(external)
        # One closed-file probe establishes the hole-punch representation before
        # generation. Every later file is checked after close, with a conservative
        # pre-write bound that includes full allocation of the next file.
        def create(path,size):
            with path.open("wb") as stream:
                stream.truncate(size)
                stream.write(b"S" * min(size,BLOCK))
                stream.flush()
                if size > BLOCK:
                    fcntl.fcntl(stream.fileno(),99,struct.pack("IIqq",0,0,BLOCK,
                        math.ceil((size-BLOCK)/BLOCK)*BLOCK))
            st=path.stat()
            assert st.st_size == size and st.st_blocks*512 <= BLOCK,(path,st.st_size,st.st_blocks)
            return st.st_blocks*512

        probe=fixture / "allocation-probe"
        create(probe,1024**2)
        probe.unlink()
        allocated=0
        peak_prewrite_bound=0
        overhead_allowance=16*1024**2
        start=time.monotonic()
        for index,(relative,size,mode) in enumerate(entries,1):
            next_bound=allocated+max(BLOCK,math.ceil(size/BLOCK)*BLOCK)+overhead_allowance
            assert next_bound <= LIMIT,"stop before allocation limit"
            peak_prewrite_bound=max(peak_prewrite_bound,next_bound)
            path=target/relative
            path.parent.mkdir(parents=True,exist_ok=True)
            allocated+=create(path,size)
            if index%50000==0:
                print(json.dumps(dict(generated=index,regular_allocated_bytes=allocated)),flush=True)
        initial_census=census(target)
        assert initial_census["regular_files"] == 393568  # native root .git pointer
        assert initial_census["directories"] == 6794
        assert initial_census["symlinks"] == 1
        entire=census(fixture)
        assert entire["allocated_regular_bytes"] <= LIMIT
        assert entire["allocated_regular_bytes"]-allocated < overhead_allowance
        result["representative"]=dict(initial=initial_census,entire_fixture=entire,
            generation_seconds=time.monotonic()-start,generated_regular_allocated_bytes=allocated,
            conservative_peak_regular_allocation_bound=peak_prewrite_bound,
            allocation_bound_satisfied=True,
            representation="Exact paths, logical sizes and directory shape; synthetic first-block bytes with punched holes; eight manifest symlinks are regular stand-ins. Additional .lake symlink targets only an owned external sentinel.")
        start_calls=len(calls)
        before=resource.getrusage(resource.RUSAGE_CHILDREN)
        start=time.monotonic()
        with patch.object(owner.subprocess,"run",side_effect=measured):
            owner.remove_registered(dummy,target)
        wall=time.monotonic()-start
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        assert wall < 300
        final_state=state(protocol,dummy,target)
        assert final_state == dict(root=False,pointer=False,registered=False)
        assert external.read_bytes()==b"owned external sentinel\n"
        assert refs(dummy)==before_refs
        result["representative"].update(wall_seconds=wall,
            user_cpu_seconds=after.ru_utime-before.ru_utime,
            system_cpu_seconds=after.ru_stime-before.ru_stime,
            native_calls=calls[start_calls:],final_state=final_state,refs_preserved=True,
            external_target_preserved=True,deadline_seconds=300)
        assert bindings(source)==initial,"source binding changed"
        result.update(passed=True,all_dummy_ref_oids=before_refs,source_revalidated=True,
            limitations=["Single candidate run under shared-host conditions; no general speed ratio or all-workload guarantee.",
            "Synthetic empty-index fixture; no tracked-index or source-content throughput claim.",
            "Sparse logical sizes preserve namespace and first-block allocation, not historical physical layout.",
            "Allocation bound covers regular st_blocks*512, not filesystem-internal metadata or shared-volume usage.",
            "Historical storage cause and registration-loss timing/actor remain unknown."])
    finally:
        start=time.monotonic()
        shutil.rmtree(fixture)
        result.update(cleanup_seconds=time.monotonic()-start,fixture_removed=not fixture.exists(),
                      owned_jobs_running=0)
        output=Path(__file__).with_name("result.json")
        output.write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps(result),flush=True)


if __name__ == "__main__":
    main()
