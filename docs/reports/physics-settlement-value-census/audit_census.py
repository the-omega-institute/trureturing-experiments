#!/usr/bin/env python3
"""Audit the fixed 36-row census against PR labels and immutable git sources.

Reads git objects, never checks out or changes the main repository. --pr-cache
may supply previously fetched `gh pr view --json` responses; these process
copies are not published. The output retains only field evidence and sources.
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import subprocess

SNAPSHOT = "53ec26edbed43dc4f4862d2db303ad80cb3dbd92"
SOURCE_URL = "https://github.com/the-omega-institute/trureturing"
HERE = Path(__file__).resolve().parent
# Reviewed locators in the PR bodies, not a semantic proof-shape classifier.
SHAPE_LINES = {
    11696: 88, 12188: 92, 12833: 88, 13018: 9, 13245: 39, 13409: 90,
    13680: 93, 14064: 230, 14159: 92, 14605: 142, 15193: 74, 15354: 146,
    15530: 233, 15565: 105, 16103: 76, 16361: 29, 16702: 66, 16709: 100,
    16727: 67, 16879: 99, 16887: 75, 16902: 63, 16949: 86, 16971: 66,
    17103: 192, 17145: 87, 17178: 53, 17487: 66, 17489: 65, 17605: 69,
    17628: 69, 17653: 93, 17760: 97, 17778: 66, 17791: 181,
}
ADMISSION_LINES = {
    11696: [92], 12188: [72], 12833: [88], 13018: [10], 13245: [35],
    13409: [78], 13680: [95], 14064: [224], 14159: [64], 14605: [142],
    15193: [74], 15354: [78], 15530: [53], 15565: [105], 15993: [35],
    16103: [70, 74], 16361: [18], 16702: [65], 16709: [102], 16727: [69],
    16879: [99], 16887: [93], 16902: [51], 16949: [88], 16971: [56],
    17103: [196], 17145: [76], 17178: [51], 17487: [58], 17489: [57],
    17605: [59], 17628: [61], 17653: [105], 17760: [101], 17778: [60],
    17791: [311, 312, 313, 314, 315],
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def without_comments(text):
    """Mask nested Lean comments and strings, preserving offsets and lines."""
    out, position, depth = [], 0, 0
    outside = re.compile(r'/-|--[^\n]*|"(?:\\.|[^"\\])*"', re.S)
    inside = re.compile(r'/-|-/')
    while position < len(text):
        match = (inside if depth else outside).search(text, position)
        end = match.start() if match else len(text)
        segment = text[position:end]
        out.append(re.sub(r'[^\n]', ' ', segment) if depth else segment)
        if not match:
            break
        token = match.group()
        if token == '/-':
            depth += 1
        elif token == '-/':
            depth -= 1
        out.append(re.sub(r'[^\n]', ' ', token))
        position = match.end()
    if depth:
        raise ValueError("unclosed source comment")
    return ''.join(out)


class GitSource:
    def __init__(self, repo):
        self.repo = str(repo)

    def git(self, *args):
        return subprocess.run(['git', '-C', self.repo, *args], capture_output=True, check=True).stdout

    def read(self, path):
        return self.git('show', SNAPSHOT + ':' + path).decode()

    def ref(self, path, text=None, line=None, excerpt=None):
        if text is None:
            text = self.read(path)
        result = dict(kind='git-file', commit=SNAPSHOT, path=path,
                      sha256=sha256(text.encode()), url=f'{SOURCE_URL}/blob/{SNAPSHOT}/{path}')
        if line is not None:
            result.update(line=line, excerpt=excerpt if excerpt is not None else text.splitlines()[line - 1])
            result['url'] += f'#L{line}'
        return result


def collect_sources(source, settlements):
    added = {}
    for row in settlements:
        merge = row['merge_commit_sha']
        added[row['pr']] = [p for p in source.git('diff', '--name-only', '--diff-filter=A',
                             '--no-renames', merge + '^1', merge, '--', 'D5').decode().splitlines()
                            if p.endswith('.lean')]
    wanted = {p for paths in added.values() for p in paths}
    entries = source.git('ls-tree', '-rz', SNAPSHOT, '--', 'D5').split(b'\0')
    entries = [e.split(b'\t', 1) for e in entries if e]
    entries = [(meta.split()[2].decode(), path.decode()) for meta, path in entries
               if path.endswith(b'.lean')]
    process = subprocess.Popen(['git', '-C', source.repo, 'cat-file', '--batch'],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    delivered_texts, edges = {}, []
    for oid, path in entries:
        process.stdin.write((oid + '\n').encode())
        process.stdin.flush()
        header = process.stdout.readline().decode().split()
        if len(header) != 3 or header[1] != 'blob':
            raise ValueError(f'cannot read source blob {oid}')
        body = process.stdout.read(int(header[2]))
        if process.stdout.read(1) != b'\n':
            raise ValueError('git batch framing')
        text = body.decode()
        if path in wanted:
            delivered_texts[path] = text
        masked = without_comments(text)
        # All import declarations in this fixed snapshot list modules on their
        # declaration line. Check that premise instead of guessing continuations.
        if re.search(r'(?m)^[ \t]*(?:(?:public|private|meta)[ \t]+)*import[ \t]*$', masked):
            raise ValueError(f'multiline import needs inspection: {path}')
        for match in re.finditer(r'(?m)^[ \t]*(?:(?:public|private|meta)[ \t]+)*import[ \t]+([^\n]*)', masked):
            for module in re.finditer(r'\bD5\.[\w.]+', match.group(1)):
                offset = match.start(1) + module.start()
                line = text.count('\n', 0, offset) + 1
                edges.append(dict(path=path, line=line, target=module.group(),
                                  excerpt=text.splitlines()[line - 1].strip(), blob=oid))
    process.stdin.close()
    if process.wait() != 0:
        raise ValueError('git source scan failed')
    return added, delivered_texts, edges, len(entries)


def checked(value, evidence, basis, status='verified'):
    return dict(value=value, status=status, basis=basis, evidence=evidence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--main-repo', required=True, type=Path)
    parser.add_argument('--pr-cache', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    labels = json.loads((HERE / 'labels.json').read_text())
    original_batch = json.loads((HERE / 'physics-settlements-36.batch.json').read_text())
    assert sha256((HERE / 'physics-settlements-36.batch.json').read_bytes()) == \
        '8a26793329161a9b0654d5976e187c8e2d167b72cebd273fab57e51d6d90004d'
    settlements = labels['settlements']
    assert len(settlements) == 36
    source = GitSource(args.main_repo)
    added, texts, edges, file_count = collect_sources(source, settlements)
    original = original_batch[0]['finite']
    roles = [x['name'] for x in original_batch[0]['plan']['roles']]
    audited, corrections = [], []
    for i, row in enumerate(settlements):
        pr = row['pr']
        if args.pr_cache:
            pr_data = json.loads((args.pr_cache / f'{pr}.json').read_text())
        else:
            pr_data = json.loads(subprocess.run(['gh', 'pr', 'view', str(pr), '--repo',
                'the-omega-institute/trureturing', '--json', 'number,title,body,url,mergedAt,mergeCommit'],
                capture_output=True, check=True).stdout)
        assert pr_data['number'] == pr and original['states'][i]['pr'] == pr
        assert pr_data['mergeCommit']['oid'] == row['merge_commit_sha']
        assert pr_data['mergedAt'][:10] == row['merged_at']
        body = pr_data['body'].splitlines()

        def pr_ref(number=None, excerpt=None, field=None):
            result = dict(kind='pr-body' if number else 'pr-metadata', url=row['url'])
            if number:
                selected_text = body[number - 1] if excerpt is None else excerpt
                result.update(body_line=number, excerpt=selected_text[:600],
                              excerpt_truncated=len(selected_text) > 600,
                              body_sha256=sha256(pr_data['body'].encode()))
            else:
                result.update(field=field, value=pr_data[field])
            return result

        prior = {k: original['outputs'][j][i]['value'] for j, k in enumerate(roles)}
        assert all(prior[k] == row[k] for k in ['shape', 'admission', 'resolution', 'track'])
        fields = {}
        if pr == 15993:
            path = added[pr][0]
            line = next(n for n, s in enumerate(texts[path].splitlines(), 1) if 'result: bind-only;' in s)
            fields['shape'] = checked('bind-only', [source.ref(path, texts[path], line)],
                'explicit source-header classification; original PR result label was inferred')
        else:
            evidence = pr_ref(SHAPE_LINES[pr])
            assert row['shape'] in evidence['excerpt'], (pr, 'shape locator changed')
            fields['shape'] = checked(row['shape'], [evidence], 'explicit PR classification, checked as a source label')
        admission_evidence = [pr_ref(n) for n in ADMISSION_LINES[pr]]
        combined = '\n'.join(x['excerpt'] for x in admission_evidence)
        for term in ['escape-witness', 'open-problem-resolution']:
            if term in row['admission']:
                assert term in combined, (pr, 'admission locator changed')
        headers, utilities, resolutions, scribe_refs = [], [], [], []
        for path in added[pr]:
            text = texts[path]
            for n, line in enumerate(text.splitlines(), 1):
                if 'admission_basis:' in line:
                    headers.append(source.ref(path, text, n))
                if re.match(r'\s*utility:', line):
                    match = re.search(r'utility:\s*(?:kind=([\w-]+)|(none))', line)
                    if match:
                        utilities.append((match.group(1) or match.group(2), source.ref(path, text, n)))
            scribe = 'Blueprint/' + path[:-5] + '.scribe.cs'
            text = source.read(scribe)
            scribe_refs.append(source.ref(scribe, text))
            for n, line in enumerate(text.splitlines(), 1):
                match = re.search(r'ResolutionKind\.(Proved|Refuted)', line)
                if match:
                    resolutions.append((match.group(1), source.ref(scribe, text, n)))
        fields['admission'] = checked(row['admission'], admission_evidence + headers,
            'PR basis for newly added modules, cross-checked against present source headers; existing owner re-pins excluded')
        kinds = {x[0] for x in utilities}
        utility = next((k for k in ['certified-instance', 'checker', 'none'] if k in kinds), None)
        fields['utility'] = checked(utility, [x[1] for x in utilities],
            'strongest explicit utility header over modules added by the merge', 'verified' if utility else 'missing')
        resolution_values = {x[0] for x in resolutions}
        assert len(resolution_values) <= 1, (pr, 'multiple conclusion types')
        resolution = next(iter(resolution_values)) if resolution_values else 'none (library result)'
        assert resolution == row['resolution'], (pr, 'conclusion mismatch', resolution)
        fields['resolution'] = checked(resolution, [x[1] for x in resolutions] or
            [pr_ref(5 if pr == 13018 else 18)] + scribe_refs,
            'explicit Scribe resolution kind' if resolutions else
            'no external resolution claim in all added-module Scribes; library scope in PR')
        track_hits = [(n, s) for n, s in enumerate(body, 1)
                      if re.search(r'\b(quick|deep)[ -]track\b', s, re.I)]
        if row['track'] == 'unlabelled':
            assert not track_hits, (pr, 'unexpected explicit track')
            fields['track'] = checked(None,
                [dict(kind='pr-body-search', url=row['url'], body_sha256=sha256(pr_data['body'].encode()),
                      pattern='explicit quick/deep track label', matches=0)],
                'no explicit quick/deep track in PR body; seat suffix is not imputed', 'missing')
        else:
            hits = [(n, s) for n, s in track_hits if re.search(r'\b' + row['track'] + r'[ -]track\b', s, re.I)]
            assert hits, (pr, 'track source missing')
            fields['track'] = checked(row['track'], [pr_ref(hits[0][0])], 'explicit PR search-track attribution')
        date = pr_data['mergedAt'][:10]
        period = 'early (2026-10-01..09)' if date <= '2026-10-09' else 'late (2026-10-10..16)'
        fields['period'] = checked(period, [pr_ref(field='mergedAt')], 'merge timestamp bucket')
        delivered = {p[:-5].replace('/', '.') for p in added[pr]}
        external_edges = [e for e in edges if e['target'] in delivered and
                          e['path'][:-5].replace('/', '.') not in delivered]
        external = sorted({e['path'][:-5].replace('/', '.') for e in external_edges})
        downstream = 'has external D5 importer' if external else 'no external D5 importer'
        fields['downstream'] = checked(downstream,
            [dict(kind='snapshot-import-scan', commit=SNAPSHOT, d5_lean_files=file_count,
                  delivered_modules=sorted(delivered), external_importers=external)] +
            [dict(kind='git-import', commit=SNAPSHOT, **e,
                  url=f"{SOURCE_URL}/blob/{SNAPSHOT}/{e['path']}#L{e['line']}") for e in external_edges],
            'external D5 source import of an added module; separate from live proof-path use')
        assert sorted(delivered) == original['states'][i]['delivered_modules'], (pr, 'delivered module mismatch')
        if external != original['states'][i]['external_importers']:
            corrections.append(dict(pr=pr, field='external_importers',
                original=original['states'][i]['external_importers'], checked=external,
                added=sorted(set(external) - set(original['states'][i]['external_importers'])),
                removed=sorted(set(original['states'][i]['external_importers']) - set(external)),
                cause='original build_census.py reads only the first 30000 characters per importer; '
                      'BroadcastingChannelExtremaRefutation import is at zero-based character 36829, line 309'
                      if pr == 12833 else 'see field source evidence'))
        for name, value in fields.items():
            if prior[name] != value['value']:
                corrections.append(dict(pr=pr, field=name, original=prior[name], checked=value['value'],
                                        status=value['status']))
        audited.append(dict(pr=pr, url=row['url'], merge_commit_sha=row['merge_commit_sha'],
            merged_at=pr_data['mergedAt'], original=prior,
            original_label_provenance=dict(shape='inferred from helpers and exact computation; no explicit result classification in PR'
                if pr == 15993 else 'PR self-report', admission='PR self-report',
                resolution='inferred from PR tightness answer' if pr == 14159 else 'PR body',
                track='missing' if row['track'] == 'unlabelled' else 'PR self-report',
                utility='source-derived', period='merge metadata', downstream='source import proxy'),
            checked=fields, delivered_modules=sorted(delivered), external_importers=external,
            live_proof_path_use=dict(status='not-audited', value=None),
            original_shape_inference_evidence=[pr_ref(39), pr_ref(40)] if pr == 15993 else []))
    counts = collections.Counter(v['status'] for r in audited for v in r['checked'].values())
    audit = dict(version=2, source_snapshot=SNAPSHOT,
        input_experiment_commit='c42809bac1b923b5fc78aaa669125516c3a65ac7',
        source_scope=dict(d5_lean_files=file_count, import_edges=len(edges),
            import_graph_sha256=sha256(json.dumps(edges, sort_keys=True).encode()),
            import_scan='full blobs, nested comments and string literals masked; source header import syntax; no length cutoff'),
        verification_boundary='verified means source-label or repository-read provenance was checked; '
            'no independent semantic proof-shape, admission, utility or Lean-kernel re-verification is claimed',
        missingness_boundary='shared missing marker means shared missingness only, not equality of unknown real search tracks',
        summary=dict(rows=36, fields_per_row=7, cells=dict(counts),
            fully_verified_rows=[r['pr'] for r in audited if all(v['status'] == 'verified' for v in r['checked'].values())],
            missing_fields=[dict(pr=r['pr'], field=k) for r in audited for k, v in r['checked'].items() if v['status'] == 'missing'],
            inferred_fields=[dict(pr=r['pr'], field=k) for r in audited for k, v in r['checked'].items() if v['status'] == 'inferred'],
            original_inferred_shape_rows=[15993],
            original_inference_now_explicitly_source_labelled=[15993]),
        corrections=corrections, settlements=audited)
    write(args.out / 'labels-audited-v2.json', audit)
    corrected = []
    for request in original_batch:
        request['binding']['label'] = 'audited v2: ' + request['binding']['label']
        request['binding']['audit'] = 'labels-audited-v2.json'
        request['binding']['boundary'] += '; missing tracks remain missing; label verification is source provenance only'
        meanings = {
            'shape': 'explicit PR proof_shape label; #15993 checked against its explicit source header',
            'track': 'explicit quick/deep PR track; null with status=missing otherwise; shared missingness is not actual track equality',
            'downstream': 'TASK: a D5 module outside the added-module set imports an added module in the fixed source snapshot; not live proof-path use',
        }
        request['binding']['readout_sources'].update(meanings)
        for role in request['plan']['roles']:
            if role['name'] in meanings:
                role['meaning'] = meanings[role['name']]
        finite = request['finite']
        finite['outputs'] = [[dict(value=r['checked'][k]['value'], status=r['checked'][k]['status'])
                              for r in audited] for k in roles]
        finite['task_outputs'] = finite['outputs'][request['plan']['task']]
        n = len(audited)
        finite['kernels'] = [[[values[i] == values[j] for j in range(n)] for i in range(n)]
                             for values in finite['outputs']]
        plan = request['plan']
        finite['layers'] = []
        for position in range(len(plan['additions']) + 1):
            selected = plan['initial'] + [k for a in plan['additions'][:position] for k in a]
            finite['layers'].append([[all(finite['kernels'][k][i][j] for k in selected)
                                      for j in range(n)] for i in range(n)])
        for state, row in zip(finite['states'], audited):
            state['external_importers'] = row['external_importers']
        corrected.append(request)
    (args.out / 'physics-settlements-36.batch.v2.json').write_text(
        json.dumps(corrected, ensure_ascii=False, separators=(',', ':')) + '\n')
    print(json.dumps(dict(summary=audit['summary'], corrections=corrections, source_scope=audit['source_scope']), ensure_ascii=False))


if __name__ == '__main__':
    main()
