# E2a prospective protocol

## Claim and authority

This protocol implements the fourth-panel two-slot pilot. It tests only process
feasibility and trustworthiness of readings for L1 autonomous-research and L2
research_readouts. It makes no decision about selection-policy effectiveness,
registers no value function f, expands no sample, and implies no long-run value.
CLAUDE.md and the specification at the pinned main baseline retain authority;
§§3.2–3.3, existing admission, and the formal three-seat review are unchanged.
This private repository holds the protocol, frame, necessary programs/data and
findings. No process logs or command streams are retained here.

## Immutable inputs and versions

Main repository: `the-omega-institute/trureturing`, commit
`442e8d0a8855794547208a39a9ae17cd9798e572`.
Read-only checkout:
`/Users/auricstudio/Library/Application Support/heca/worktrees/568cd92c/lively-pebble`.
Experiment base: `a01eb2804313522286824a4787052c4ab630be11`.
Branch: `exp/autonomous-research-e2a-1017`.

The following Git blob identities bind versions, not cache eligibility:

| Component | Git blob |
| --- | --- |
| skills/autonomous-research/SKILL.md (L1) | bd4ea908ba8ed2f3153a7895a4c3d126b9b48f7b |
| tools/scripts/agent/research_readouts/readouts.py (L2) | 5e4d264cb223a2c8287e5472179fccfaca887c09 |
| tools/scripts/agent/research_readouts/README.md | cf50fae2514789482383a7f73df1dafcc833005a |
| tools/scripts/agent/research_readouts/test_readouts.py | b745f4e5680d20d5934050f2a820aa54680056c2 |
| tools/lean-inspector/LeanInformationAudit/CompiledAxioms.lean | b023885850622eb8708c4332391ade252d451c0c |
| CLAUDE.md | bfd12658202af15ac014c0e765b86efb7845ec85 |
| docs/develop/spec/golden-ledger-repo-spec.md | 76df18f3fd062092c545288cc80fe020f0d00c3e |
| lean-toolchain | 025e59548e48cf71f2744154e2890beefe30a258 |
| lake-manifest.json | 7f95ee71abeba821317bb90812b504097b2abc7e |

Enumeration uses Python standard library (observed Python 3.9.6) and Git
(observed 2.50.1), with source bytes read through the immutable Git tree.
No date, working-directory path, clock reading or runtime duration enters
canonical frame output. The acquisition record separately describes available
compiled inputs. The frame program and its first outputs follow this protocol
in Git history. This document must be committed before running either frame.

## Generation budgets and non-probe boundary

Two sources, in order: boundary, then residual. No substitution or third source.
Boundary generation: at most 600 wall seconds from its first enumerator
invocation, one immutable source inventory, one capability inspection of the
pinned checkout's existing compiled inputs. Residual generation: at most 1800
wall seconds from its first enumerator invocation, one complete immutable source
inventory, and at most the first 20 registrations in sorted order inspected for
an exact event and known-result status. Byte reproducibility reruns of the same
inventory are allowed; they cannot change the candidate rules or inputs.
Each inspected prospective event has at most 120 seconds of repository novelty
lookup and 180 seconds of literature metadata/statement lookup, within the
source's wall budget. A network failure is a missing reading, not no hit.
An entry left after the inspected prefix is explicitly not assessed, not skipped.
Stop on the first eligible event. A missing prerequisite, exhausted budget or
no eligible event makes that source infeasible and consumes its original slot.
The other original source may still be described; it cannot repair the failed
slot or permit an affirmative two-slot feasibility claim.

The proposer does not solve: no proof search, counterexample search, tactic
execution, model enumeration, arithmetic trials, symbolic answer computation or
mathematical sampling. Dependency extraction and reading statements/metadata
are allowed only to construct the frame. Known-result screening compares
already-stated results, their assumptions and scope; it never requests or
constructs a proof/refutation. Do not follow a source into a proof. If metadata
reveals the answer, exclude the duplicate before freezing, record that source,
and continue only within the fixed original generation budget. Any mathematical
probe by the proposer invalidates the affected slot and is reported; no repair
by replacement is allowed.

## Boundary frame

The universe is the public theorem declarations of D5 modules frozen at the
immutable baseline, as identified by Golden/Frozen/state and current baseline
source. For each theorem obtain distinct external declaration users on actual
compiled proof-value dependency paths. External means outside the theorem's
defining module; generated implementation constants may carry paths but do not
inflate the distinct user count. List the named users and their paths. Type-only
references, imports, source-name hits and header consumer assertions are not
proof-use evidence.

CompiledAxioms.declarationDependencyParts separates type/constructor references
from optional proof/definition value references. A value-edge path must start
at the external declaration's actual value and reach the theorem; type edges
are excluded. This is a compiled syntactic dependency reading, not a claim
about semantic indispensability or §3.2 proof shape. The compiled environment
must cover the complete baseline D5 user domain before claiming exact counts
and sorting; partial coverage is missing, not zero. Lean's elaborated constants,
module ownership and private/generated-name classification form the acquisition
trust boundary. Kernel validity does not itself certify the extraction procedure.

For every theorem with at least two verified distinct external users, generate
one candidate per explicit Prop hypothesis removable on the original carrier.
Retain all other parameters, hypotheses and the conclusion, with their original
quantifiers. Removal must leave a well-formed telescope and conclusion; no
remaining binder/type may depend on the removed proof variable. Do not alter
typeclasses or carriers. Sort by verified user count descending, fully qualified
theorem name ascending, original telescope hypothesis index ascending (one-based,
counting all binders). Exclude already-known targets using the fixed duplicate
screen, and freeze the first eligible row. Do not let a solver choose.

If complete declaration-level proof-use readings cannot be acquired without
writing/building the read-only checkout, stop this source as infeasible. Report
its frozen source inventory, unknown counts and unavailable candidate expansion;
do not rank unknowns as zero, use import counts, scrape textual proof references
as a proxy, build elsewhere, or obtain a substitute population.

## Residual frame

Completely enumerate registration declarations in baseline Reg/ and D5/ whose
continuation field is explicitly `.unknown`, ignoring comment/string occurrences
and distinguishing omitted/default fields. Preserve every occurrence, full name,
source line and relevant registration fields. Report unsupported syntax or an
unresolvable declaration identity as missing input; do not silently drop it.
Sort by fully qualified registration declaration name ascending, then source
path and line only to expose identity collisions (collisions must be resolved
before eligibility). Source inventory is complete even if event assessment is
budget-limited.

For each item in that order, inspect its registration, associated template,
arena, readouts and named target. Write a closed binary proposition that asks
whether those readouts determine that named target. Specify every quantifier,
parameter range, allowable operation, original arena/carrier, and common actual
implementation constraints. Do not ask merely whether residual remains, change
the arena, or independently combine incompatible realizations. `unknown` is not
evidence of nonempty residual, undecidability, novelty or any truth value.

A row can be skipped only for a source-grounded non-solving fact: a duplicate
already-stated result, an event that cannot be formed, or a missing necessary
input. Each skip records its precise source and reason. The first item with a
checkable exact event and completed duplicate screen is frozen; do not inspect
later items for a more attractive mathematical answer. If an event is not
formable within the fixed generation budget, retain the failed-source result.

## Duplicate screening and freezing

Before freezing, retain the exact repository search terms, searched pinned
paths, metadata/statement hits and bibliographic URLs or identifiers. Search
D5 including private declarations, Problems, Library, relevant theory/template
sources, and the pinned Mathlib source only if already available read-only.
Search literature indexes/publisher metadata using the source's own names and
exact relation; inspect statements only. Query open PR/issue titles and bodies
for the candidate names if access is available; keep acquisition time and scope
separate from pinned source facts. A missing Mathlib checkout or network query
is recorded as incomplete novelty evidence; do not claim exhaustive novelty.
An exact already-stated result, positive or negative, is a duplicate and is
excluded before freezing. Absence of a hit is limited to the stated search scope.
A source with unresolved identity, missing mandatory event inputs or incomplete
mandatory duplicate screening has no eligible frozen event.

Freezing retains a machine-readable statement and quantifiers, full source
identity, and exact duplicate readings. No replacement after freezing: probe
solution, all bind-only, agreement, failure to start, timeout, invalid or no target
still occupies the original slot. No solver may change the frozen text.

## Subsequent prospective execution (caller-owned)

Caller owns orchestration, final freeze and publication. Proposer owns only this
protocol, deterministic frames and events, never a mathematical probe.
Two mutually hidden prediction seats: one codex-cli and one nyxid-oracle/GPT Pro;
record actual seat/model-family availability, inherited priors and any missing
oracle response. They receive only the frozen event and allowed source context,
not each other's probabilities or solver output. Solver is a separate codex
loading the pinned autonomous-research skill, sees no probabilities, and cannot
change the target. A fourth, independent codex evaluator checks fidelity,
evidence, execution/settlement times and §3.2 shape before reading score output.
The proposer is not a predictor, solver or evaluator.

Each slot has its own main-repository preregistration issue, created later by the
caller. Predictions precede every mathematical probe and are published by caller
in the issue's unique, never-edited `research-prereg-v1` comment. This step creates
no main-repository issue, comment, branch or PR. An infeasible frame is reported
here and has no invented research event or fabricated L2 input.

For each formed unit i, let C_i be the carrier comment's GitHub REST created_at
in UTC seconds. Caller computes D_i = C_i + 3600 seconds, emits D_i as an absolute
`YYYY-MM-DDTHH:MM:SSZ` deadline in the same initial comment and verifies that its
actual created_at yields that exact deadline. Since created_at is not known
before POST, use the caller's preselected UTC publication second P_i and set
D_i = P_i + 3600 seconds in the payload; if C_i differs from P_i, the comment's
already-published absolute D_i is authoritative, and the available budget is
D_i - C_i (must be positive and at most 3600 seconds). Never edit/repost to reset
the clock. A mismatch that would give over 3600 seconds or a nonpositive window
invalidates that original slot. This rule fixes a concrete deadline before
execution without trusting a later payload edit.

Unit 1 is boundary; unit 2 is residual. Do not publish/start unit 2 before unit 1
has reached its recorded terminal, including independent evaluation. Each solver
has at most 1800 wall seconds including an initial at-most-300-second cheap
probe. Evaluator has at most 600 wall seconds, all within D_i; registration and
handoff use the remainder. At any resource/deadline trigger settle immediately,
with unresolved/invalid when applicable; do not increase budgets. A missing unit
prevents the full two-slot execution claim. Caller may document the other source
but must not relabel a one-slot observation as the planned two-slot pilot.

For each unit, publish its unchanged matching `research-terminal-v1` block with
traceable existing evidence. Only independently verified proved/refuted before
D_i supply a in {0,1}. Unresolved is never scored false. No invocation in this
proposer task supplies probabilities, answer evidence or §3.2 shape.

## Readings, success, falsification and stop

L2 scoring rule: `{"name":"binary-brier","version":1}`.
`surprise_threshold = 0.5625 = 9/16`; y=1 iff the equally weighted mean of seat
Brier losses is at least 9/16. This convention is not a value threshold and
plays no role in success. Keep missingness, denominators and actual families.
No f, E[f], quota adjustment or automatic additional unit is computed.

Process feasible requires both qualified frozen events, pre-probe predictions,
a traceable terminal by each absolute cutoff (legal unresolved counts), separated
roles and no substitution. Readings trustworthy requires canonical L2 replay
exit 0 plus evaluator's independent statement/evidence/time/shape and arithmetic
checks. At least one timely proved/refuted unit is necessary to say that the
answer-scoring chain was observed; otherwise that segment is unverified.
Selection policy is descriptive and never adjudicated by this pilot.

Inability of a source to form an event in budget, post-hoc predictions, event
change, missing terminal, role mixing or resource overrun falsifies its respective
claim. Stop/settle when any registered trigger applies. Do not expand the frame,
relax the criteria, increase a budget or manufacture a target to fill a slot.

## This delivery

Generation, checking and publication must stop before the worker's absolute
cutoff `2026-10-17T08:41:28Z`. Commit and push whatever eligible findings exist
before then; list unfinished requirements as unresolved. Only the designated
experiment worktree is writable, except the explicitly required runner result
and sentinel. Never touch `/Users/auricstudio/trureturing-experiments`.
No main-repository builds, writes, branches, issues or PRs. No force push or merge.
Open an experiment-repository PR to main, retaining the separate protocol
ancestor commit and frame/event outcome commit(s).
