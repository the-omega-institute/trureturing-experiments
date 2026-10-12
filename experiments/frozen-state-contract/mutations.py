"""Verify named contract regressions against isolated production mutations."""
import argparse
import json
from pathlib import Path
import subprocess

CLI = 'tools/tests/StrataLint.Tests/StrataLint.Tests.csproj'
FILEMAP = 'tools/tests/StrataLint.FileMap.Tests/StrataLint.FileMap.Tests.csproj'
MUTATIONS = [
    ('schema', 'tools/StrataLint.Engine/FrozenState/FrozenStateRecord.cs',
     '            RequireKeys(root, ["statement_id", "declarations"], "record");',
     '            if (root.EnumerateObject().Count() == 1) return new(StatementId.Create(root.GetProperty("statement_id").GetString()!), []);\n            RequireKeys(root, ["statement_id", "declarations"], "record");',
     CLI, 'StrictStateSchemaRejectsSingleFieldRecord'),
    ('pin', 'tools/StrataLint.Engine/Rules/RepositoryRules.FrozenState.cs',
     '            if (record.StatementId != actual)',
     '            if (false && record.StatementId != actual)',
     CLI, 'CurrentStateChecksPinAndDeclarationInventoryIndependently'),
    ('inventory', 'tools/StrataLint.Engine/Rules/RepositoryRules.FrozenState.cs',
     '            if (!record.Declarations.SequenceEqual(FrozenStateRecord.Project(report)))',
     '            if (false && !record.Declarations.SequenceEqual(FrozenStateRecord.Project(report)))',
     CLI, 'CurrentStateChecksPinAndDeclarationInventoryIndependently'),
    ('write-set', 'tools/StrataLint.Cli/Commands/Ledger/DagLedgerAlignWriter.cs',
     '            return new(true, Summary(selected.Length, changed, added, selected.Length - changed - added)',
     '            _ = FrozenStateWriter.Write(root, RepoPath.CreateKnown("D5/S0/Carrier/B.lean"), StatementId.Create("sha256:" + new string(\'f\', 64)), []);\n            return new(true, Summary(selected.Length, changed, added, selected.Length - changed - added)',
     CLI, 'RepinWritesOnlyTargetWithFrozenDownstreamAndIdenticalRetryWritesNothing'),
    ('unclassified-accepted', 'Meta/FILEMAP.toml',
     '  { pattern = "Golden/Frozen/state/**/*.json",',
     '  { pattern = "Golden/Frozen/accepted/**", kind = "data", admission_plane = "content", produced_by = "none", consumed_by = ["reader"], verified_by = ["repository-policy"], artifact_id = "none", runtime_disposition = "committed-source" },\n  { pattern = "Golden/Frozen/state/**/*.json",',
     FILEMAP, 'AcceptedPathIsRejectedByCanonicalFileMapCoverage'),
    ('release', 'tools/StrataLint.Cli/Commands/TruthReleaseCommand.cs',
     '        return new(2, string.Empty, "TRUTH_RELEASE_EVENT_FIELDS_UNAVAILABLE "',
     '        return new(0, string.Empty, "TRUTH_RELEASE_EVENT_FIELDS_UNAVAILABLE "',
     CLI, 'TruthReleaseRejectsIssuanceWhenEventFieldsHaveNoSource'),
]

def main():
    parser = argparse.ArgumentParser(); parser.add_argument('repository', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--only', nargs='+'); args = parser.parse_args()
    registrations = [{'name': name, 'mutation_position': path, 'named_red_test': test,
                      'expected_red': 'The named regression fails with zero compiler errors.'}
                     for name, path, _, _, _, test in MUTATIONS]
    expected = args.output.with_name('mutation-expectations.json')
    expected.write_text(json.dumps(registrations, indent=2)+'\n')
    results = json.loads(args.output.read_text()) if args.only and args.output.exists() else []
    for name, relative, original_text, mutated_text, project, test in MUTATIONS:
        if args.only and name not in args.only: continue
        path = args.repository / relative; original = path.read_bytes()
        text = original.decode(); assert text.count(original_text) == 1, name
        command = ['make', '-C', 'tools', 'test', 'TEST_PROJECT='+project,
                   'TEST_FILTER=FullyQualifiedName~'+test]
        try:
            path.write_text(text.replace(original_text, mutated_text))
            red = subprocess.run(command, cwd=args.repository, capture_output=True, text=True)
        finally:
            path.write_bytes(original)
        log = red.stdout+red.stderr
        compile_errors = 0 if ('0 个错误' in log or '0 Error(s)' in log) and 'error CS' not in log else None
        named_failure = test in log and '[FAIL]' in log
        green = subprocess.run(command, cwd=args.repository, capture_output=True, text=True)
        result = {'mutation_position': relative, 'name': name, 'named_red_test': test,
                  'command': command, 'compile_errors': compile_errors, 'red_exit': red.returncode,
                  'named_failure': named_failure, 'restored': path.read_bytes() == original,
                  'green_exit': green.returncode, 'preregistered_expected_red': next(item['expected_red'] for item in registrations if item['name'] == name)}
        results = [item for item in results if item['name'] != name]
        results.append(result); args.output.write_text(json.dumps(results, indent=2)+'\n')
        print(json.dumps(result), flush=True)
        if red.returncode == 0 or compile_errors != 0 or not named_failure or green.returncode != 0:
            raise SystemExit('mutation evidence incomplete: '+name)

if __name__ == '__main__': main()
