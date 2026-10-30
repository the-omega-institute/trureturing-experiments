# Wang–Dong Conjecture 6: bounded engineering experiments

The full literal binary and arbitrary-alphabet source event is unresolved.
These are two closed bounded engineering experiments, with no surviving
counterexample, kernel proof or KPI credit. See SOURCE-LICENSE.md for scope
and rights. Fixed n=19 witness bounds constrain search, not the source event.

* `readiness/`: unrestricted table encoding, exact image oracle, initial SMT
  relaxation and measured SAT/rejection/refinement-timeout result.
* `upfront-cuts/`: the same target with both necessary one-letter cuts upfront;
  the sole check timed out without a model.
* `audit.py`, `audit-results.json`: closed independent fixed-table regression
  program and previously reported mathematical test results. Retention did
  not execute the audit or any synthesis.

Use Python 3.13 with -B and an external isolated environment containing the
pinned requirement. For the fixed-table regression, execute `audit.py` in a
disposable copy of this unit; it writes audit-results.json in the working
directory. Each experiment README describes its own reproduction entry.
Programs that perform synthesis require a separate bounded authorization.

Retention moved only the listed closed material. Original program and SMT
bytes are unchanged except audit.py's absolute input root was replaced by a
path relative to this directory. JSON projection removes UTC process anchors
and a joined-process receipt only; all mathematical and resource values remain
unchanged. Readmes and license notes were rewritten to retain substantive
scope, limitations, resources and reproduction instructions. The pilot process
ownership marker is omitted. No environments, binaries, caches or logs belong
to this unit. retention-manifest.json records hashes and these projections.
