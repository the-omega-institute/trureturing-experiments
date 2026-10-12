import json
from pathlib import Path
import tempfile
import unittest
import migrate

PIN = 'sha256:' + 'a' * 64
MODULE = 'D5/S0/Carrier/Probe.lean'
DECL = {'declaration_name_key': 'ns(n0,5:probe)', 'kind': 'theorem', 'statement_id': PIN}

class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.state = self.root / ('Golden/Frozen/state/' + MODULE + '.json')
        self.state.parent.mkdir(parents=True)
        self.state.write_text(json.dumps({'statement_id': PIN}) + '\n')
        self.accepted = self.root / 'Golden/Frozen/accepted'
        self.accepted.mkdir(parents=True)
        self.write_event('a.json')

    def write_event(self, name, pin=PIN, declarations=None):
        event = {'event_hash': PIN, 'event_type': 'Freeze', 'schema_version': 5,
                 'payload': {'descriptor_selector': MODULE, 'statement_id': pin,
                             'declaration_statement_ids': [DECL] if declarations is None else declarations,
                             'prerequisite_frozen_node_ids': []}}
        (self.accepted / name).write_text(json.dumps(event) + '\n')

    def image(self):
        return {str(p.relative_to(self.root)): p.read_bytes()
                for p in self.root.rglob('*') if p.is_file()}

    def test_projects_only_name_and_kind_and_preserves_members_and_pins(self):
        before = self.image()
        records, result = migrate.plan(self.root)
        self.assertEqual(before, self.image())
        self.assertEqual({'statement_id': PIN, 'declarations':
                          [{'declaration_name_key': DECL['declaration_name_key'], 'kind': 'theorem'}]},
                         json.loads(records[self.state]))
        self.assertEqual(1, result['members_before'])
        self.assertEqual(1, result['members_after'])
        self.assertEqual(1, result['pins_equal'])
        self.assertEqual(1, result['declaration_entries'])
        self.assertEqual(0, result['inconsistencies'])

    def test_mismatch_stops_before_writes(self):
        self.write_event('a.json', pin='sha256:' + 'b' * 64)
        before = self.image()
        with self.assertRaises(migrate.MigrationError):
            migrate.plan(self.root)
        self.assertEqual(before, self.image())

    def write_event_with_prerequisites(self, name, prerequisites, declarations=None):
        self.write_event(name, declarations=declarations)
        path = self.accepted / name
        event = json.loads(path.read_text())
        event['payload']['prerequisite_frozen_node_ids'] = prerequisites
        path.write_text(json.dumps(event) + '\n')

    def test_duplicate_freezes_with_identical_projection_collapse(self):
        self.write_event_with_prerequisites('b.json', ['sha256:' + 'c' * 64])
        records, result = migrate.plan(self.root)
        self.assertEqual({'statement_id': PIN, 'declarations':
                          [{'declaration_name_key': DECL['declaration_name_key'], 'kind': 'theorem'}]},
                         json.loads(records[self.state]))
        self.assertEqual(1, result['identical_freezes_collapsed'])
        self.assertEqual(1, result['members_after'])

    def test_duplicate_freezes_with_different_declarations_stop(self):
        other = dict(DECL, declaration_name_key='ns(n0,5:other)')
        self.write_event_with_prerequisites('b.json', [], declarations=[DECL, other])
        before = self.image()
        with self.assertRaises(migrate.MigrationError):
            migrate.plan(self.root)
        self.assertEqual(before, self.image())

    def test_duplicate_freeze_with_another_pin_stops(self):
        self.write_event('b.json', pin='sha256:' + 'b' * 64)
        with self.assertRaises(migrate.MigrationError):
            migrate.plan(self.root)

    def test_missing_freeze_stops(self):
        (self.accepted / 'a.json').unlink()
        with self.assertRaises(migrate.MigrationError):
            migrate.plan(self.root)

    def test_duplicate_declaration_stops(self):
        self.write_event('a.json', declarations=[DECL, DECL])
        with self.assertRaises(migrate.MigrationError):
            migrate.plan(self.root)

    def test_ordinal_sort_and_one_line_lf(self):
        second = dict(DECL, declaration_name_key='ns(n0,1:z)')
        self.write_event('a.json', declarations=[second, DECL])
        records, _ = migrate.plan(self.root)
        wire = records[self.state]
        self.assertEqual(1, wire.count(b'\n'))
        self.assertTrue(wire.endswith(b'\n'))
        self.assertEqual(sorted([DECL['declaration_name_key'], second['declaration_name_key']]),
                         [d['declaration_name_key'] for d in json.loads(wire)['declarations']])

if __name__ == '__main__':
    unittest.main()
