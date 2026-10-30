# Wang–Dong Conjecture 6: upfront one-letter cuts

The sole fresh solver check returned unknown (timeout); candidate is null.
The full source event remains unresolved. The fixed binary target is n=19,
m=6, s=13, r=1, S={0,...,12}, A={0}, B=17, reaching bound 126 and reset
bound 220. Both no-success conjuncts, A subset Sw AND |Sw|=13, hold as
constraints for epsilon, [0] and [1]. The formula contains 766 declared
variables and 431 assertions. These three cuts are a relaxation of the full
all-word condition through 17.

`outcome.json` preserves solver statistics and formula counts; `resources.json`
preserves measured resource data. Build took 1.3853620829759166 seconds,
serialization 1.5524347919854335 seconds and solving 60.45412520796526 seconds.
The configured timeout was 60000 ms and Z3 memory limit 1536 MB. The externally
bounded subprocess took 63.78374666697346 seconds and peaked at 960561152 bytes
RSS. The 90-second external limit did not trigger. No SAT model was available;
model decoding, admissibility replay and full oracle/certificate replay were
not exercised on n19. Their retained implementations receive no new validation
claim here.

`upfront.smt2.gz` decompresses to 2739652 bytes with SHA256
`a9948c015d90ab55cceb8223a2e5707af26dcf1f9d2d977d62437a0c0293aa11`.
`encoding.py`, `experiment.py` and `supervisor.py` retain the one-check program.
For reproduction use Python 3.13, install `requirements.txt` in an environment
outside this repository, then run `python3.13 -B supervisor.py` from a disposable
copy of this directory. This invokes one 60-second solver check under the
90-second process bound and overwrites output files in that copy.

Unknown, bounded UNSAT, finite nonfindings and restricted witness bounds have
no universal proof force. A future survivor needs independent replay and
literal kernel/admission closure before any full-resolution claim.
