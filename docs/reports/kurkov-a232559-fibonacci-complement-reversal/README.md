# Kurkov’s A232559 Fibonacci complement-reversal conjecture

Run `python3 check.py` with Python 3.10 or later. The program uses only the standard library and exact integer arithmetic; it needs no network access. `check.out` contains the output.

The fidelity checks reproduce every displayed DATA term in OEIS A232559 revision 63 (66 terms), A059894 revision 97 (71), A348366 revision 23 (65), A072649 revision 104 (99), and A130312 revision 39 (88), as well as the first five generations printed in A232559 COMMENTS. The embedded lists are transcriptions of those published sequence data, attributed to the OEIS Foundation and the entry contributors.

The identity check compares a literal FIFO generator with immediate duplicate removal against the independently implemented A348366 recurrence followed by binary complement-reversal, for all indices 1 through 50,000. A separate enumeration of binary suffix words of weighted cost at most 21, sorted by their operation paths, checks 46,367 nodes: source order, Fibonacci row lengths and index intervals, prefix decomposition, and complement-reversal conjugation. These are finite sanity checks and do not establish the all-index conjecture or any Lean proof.

Sources: [A232559](https://oeis.org/A232559), [A059894](https://oeis.org/A059894), [A348366](https://oeis.org/A348366), [A072649](https://oeis.org/A072649), [A130312](https://oeis.org/A130312). OEIS data are available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The original checking code is released under the MIT license.
