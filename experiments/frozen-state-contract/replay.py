"""Replay historical merges using canonical module state records and source imports."""
import argparse
from collections import defaultdict
import json
import hashlib
import os
from pathlib import Path
import re
import subprocess
import tempfile

from migrate import canonical, project

CASES = ['627a558a9c', '09674ad841', 'f8137cae8d', '9a29946a31',
         'c1321c1f62', '146d4f8678', '57bab5c18b', '1f0c36bb2c',
         '756b195ca3', 'e2f35ca445', '7284119321', '6d3dcd868c', '2b7dbda6e5']
STATE = 'Golden/Frozen/state/'
ACCEPTED = 'Golden/Frozen/accepted/'

class Replay:
    def __init__(self, source, simulation):
        self.source, self.sim = str(source), str(simulation)
        self.blobs, self.trees = {}, {}
        self.env = dict(os.environ, GIT_AUTHOR_NAME='contract replay',
                        GIT_AUTHOR_EMAIL='replay@example.invalid',
                        GIT_COMMITTER_NAME='contract replay',
                        GIT_COMMITTER_EMAIL='replay@example.invalid')
        Path(self.sim).mkdir(parents=True)
        self.git(self.sim, 'init', '-q')
        objects = self.git(self.source, 'rev-parse', '--git-path', 'objects').stdout.decode().strip()
        if not os.path.isabs(objects):
            objects = str(Path(self.source) / objects)
        Path(self.sim, '.git/objects/info/alternates').write_text(objects + '\n')

    def git(self, repo, *args, data=None, allow=(0,)):
        result = subprocess.run(['git', '-C', str(repo), *args], input=data,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self.env)
        if result.returncode not in allow:
            raise RuntimeError(result.stderr.decode())
        return result

    def entries(self, revision, repo=None):
        raw = self.git(repo or self.source, 'ls-tree', '-rz', revision).stdout
        return {path.decode(): metadata.split()[2].decode()
                for row in raw.split(b'\0') if row
                for metadata, path in [row.split(b'\t', 1)]}

    def load(self, ids, repo=None):
        missing = sorted(set(ids) - self.blobs.keys())
        if missing:
            raw = self.git(repo or self.source, 'cat-file', '--batch', data=('\n'.join(missing)+'\n').encode()).stdout
            offset = 0
            for oid in missing:
                end = raw.index(b'\n', offset)
                size = int(raw[offset:end].split()[2])
                self.blobs[oid] = raw[end+1:end+1+size]
                offset = end+size+2

    def synthesize(self, revision, parent=None):
        entries = self.entries(revision)
        relevant = {path: oid for path, oid in entries.items()
                    if path.startswith((STATE, ACCEPTED))}
        self.load(relevant.values())
        freezes = defaultdict(list)
        for path, oid in relevant.items():
            if path.startswith(ACCEPTED):
                event = json.loads(self.blobs[oid])
                if event['event_type'] != 'Freeze':
                    raise ValueError('non-Freeze replay input: ' + path)
                payload = event['payload']
                freezes[payload['descriptor_selector']].append(payload)
        records, collapsed = {}, 0
        for path, oid in relevant.items():
            if not path.startswith(STATE):
                continue
            module = path[len(STATE):-5]
            pin = json.loads(self.blobs[oid])['statement_id']
            matching = [canonical(project(payload)) for payload in freezes[module]
                        if payload['statement_id'] == pin]
            if not matching or len(set(matching)) != 1:
                raise ValueError('no unique state-pin projection: ' + revision + ':' + module)
            records[path] = matching[0]
            collapsed += len(matching)-1
        self.git(self.sim, 'read-tree', revision)
        deletes = b''.join(('0 ' + '0'*40 + '\t' + path + '\n').encode() for path in relevant)
        self.git(self.sim, 'update-index', '--index-info', data=deletes)
        adds, fresh = [], {}
        for path, record in records.items():
            oid = hashlib.sha1(b'blob '+str(len(record)).encode()+b'\0'+record).hexdigest()
            if oid not in self.blobs:
                fresh[oid] = record
            self.blobs[oid] = record
            adds.append('100644 '+oid+'\t'+path+'\n')
        with tempfile.TemporaryDirectory() as scratch:
            paths = []
            for oid, record in fresh.items():
                file = Path(scratch, oid); file.write_bytes(record); paths.append(str(file))
            if paths:
                self.git(self.sim, 'hash-object', '-w', '--stdin-paths', data=('\n'.join(paths)+'\n').encode())
        self.git(self.sim, 'update-index', '--index-info', data=''.join(adds).encode())
        tree = self.git(self.sim, 'write-tree').stdout.decode().strip()
        commit = self.git(self.sim, 'commit-tree', tree, '-m', revision,
                          *(['-p', parent] if parent else [])).stdout.decode().strip()
        return commit, {'members': len(records), 'identical_freezes_collapsed': collapsed}

    def verify_tree(self, tree):
        entries = self.entries(tree, self.sim)
        states = {path: oid for path, oid in entries.items() if path.startswith(STATE)}
        sources = {path: oid for path, oid in entries.items() if path.endswith('.lean')}
        # Read merged blobs from the simulation object store, including conflict-marker records.
        missing, imports = [], 0
        module_names = {path[:-5].replace('/', '.') for path in sources}
        for path in states:
            module = path[len(STATE):-5]
            if module not in sources:
                missing.append({'module': module, 'reason': 'source absent'})
        self.load(sources.values(), self.sim)
        for path, oid in sources.items():
            if not path.startswith(('D5/', 'Reg/')):
                continue
            text = self.blobs[oid].decode()
            text = re.sub(r'/\-.*?\-/', '', text, flags=re.S)
            for line in text.splitlines():
                line = line.split('--', 1)[0].strip()
                if line.startswith('import '):
                    for name in line[7:].split():
                        if name.startswith(('D5.', 'Reg.', 'Interface.', 'Impl.')):
                            imports += 1
                            if name not in module_names:
                                missing.append({'module': path, 'import': name})
        return {'records_unique_by_path': len(states) == len(set(states)),
                'record_paths': len(states), 'source_import_edges_checked': imports,
                'missing_source_imports': missing,
                'accepted_paths': sum(path.startswith(ACCEPTED) for path in entries)}

    def run(self, case):
        ours, theirs = case+'^1', case+'^2'
        base = self.git(self.source, 'merge-base', ours, theirs).stdout.decode().strip()
        b, bc = self.synthesize(base)
        o, oc = self.synthesize(ours, b)
        t, tc = self.synthesize(theirs, b)
        merged = self.git(self.sim, '-c', 'merge.renames=false', 'merge-tree',
                          '--write-tree', '--name-only', o, t, allow=(0, 1))
        lines = merged.stdout.decode().split('\n\n', 1)[0].splitlines()
        return {'case': case, 'ours': ours, 'theirs': theirs, 'base': bc,
                'ours_projection': oc, 'theirs_projection': tc,
                'merge_exit': merged.returncode, 'conflict_paths': lines[1:],
                'merged_tree': self.verify_tree(lines[0])}

def direct_merge(root, same):
    env = dict(os.environ, GIT_AUTHOR_NAME='replay', GIT_AUTHOR_EMAIL='r@example.invalid',
               GIT_COMMITTER_NAME='replay', GIT_COMMITTER_EMAIL='r@example.invalid')
    def git(*args):
        return subprocess.run(['git','-C',str(root),*args], capture_output=True, env=env)
    root.mkdir(parents=True)
    git('init','-q')
    path = root / (STATE+'D5/S0/Carrier/Probe.lean.json')
    path.parent.mkdir(parents=True)
    def commit(pin):
        path.write_bytes(canonical({'statement_id':'sha256:'+pin*64, 'declarations':[]}))
        git('add','-A'); result=git('commit','-qm','state'); assert result.returncode==0
        return git('rev-parse','HEAD').stdout.decode().strip()
    base=commit('a'); ours=commit('b'); git('checkout','-q','--detach',base)
    theirs=commit('b' if same else 'c')
    merged=git('-c','merge.renames=false','merge-tree','--write-tree','--name-only',ours,theirs)
    return {'case':'same-content' if same else 'different-records', 'merge_exit':merged.returncode,
            'conflict_paths':merged.stdout.decode().split('\n\n',1)[0].splitlines()[1:]}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('repository',type=Path)
    parser.add_argument('--output',type=Path,required=True); args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='frozen-state-replay-') as temporary:
        root=Path(temporary)
        replay=Replay(args.repository,root/'history')
        results=[replay.run(case) for case in CASES]
        direct=[direct_merge(root/'different',False),direct_merge(root/'same',True)]
    document={'historical':results,'direct':direct}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(document,indent=2)+'\n')
    print(json.dumps({'historical_cases':len(results),'historical_conflicts':sum(r['merge_exit']==1 for r in results),
                      'missing_source_imports':sum(len(r['merged_tree']['missing_source_imports']) for r in results),
                      'direct':direct}))

if __name__ == '__main__': main()
