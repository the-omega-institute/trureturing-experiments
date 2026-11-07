# Bounded-rank IID cutoff checks

Source: Lionel J. Dmello, Xiangling Xu, Marc-Olivier Renou and David Gross,
*An operational characterization of finite-dimensional quantum theory*,
arXiv:2610.03694v2, Sections II–III and VII.

Run `python3 check.py` with Python 3 and SymPy. All evidence uses integer,
rational or exact algebraic arithmetic. The expected final line is `ALL_CHECKS_OK`.
The captured output is `check.out`.

The program derives equation (1) from two Bell states and Erwin's rotated
Bell-basis effects for all twelve settings/outcomes at cos(theta)=4/5,
sin(theta)=3/5. It verifies the sign matrices' independence, measurement sum,
chronological order, normalized targets and Lemma 5 probabilities on every
zero-padded word of length at most four, at two exact rotations (5,690 words).
It also checks 96 uniquely solvable finite normalization steps.

For the proposed D-3 implication it scans 36 explicit linear datasets with
rank D=4,...,9, perturbation parameter -1,0,1, and two hidden-dynamics forms,
using exact reachable spans through length 2D. Every dataset passing the
cutoff passes the scan. These checks support the conventions and the proof
outline; they do not prove the universally quantified theorem.

Triage: the delay-line datasets first fail at D-3 while passing through D-4.
They obey the measurement-completeness identity and have fixed maps across
positions and lengths, but their realization by GPT state/effect cones has
not been established. Thus they show sharpness only in the relaxed linear
model. Sharpness among physical GPT realizations remains open.

The planned positive settlement uses a general two-sided span and annihilator
argument. Finite enumeration is not its main mathematical content.
