# Independent monophonic product bound

Run the exact check with Python 3.12 or later.

    python3 check.py

The program enumerates colorings and induced paths for complete and edgeless graphs of orders 1 through 6, reproducing instances of Lemma 3.5. It then checks induced five-vertex paths for all triples of word vertices, explicit pair colorings, and connectivity on both sides of complementation for the binary-incidence graph at k = 3, 4, 5, 6. The lower bound gives strict violations at k = 5 and k = 6. A separate calculation checks the product margin and recurrence through k = 100. The finite readings are in check.out. [Source, arXiv:2408.13494v2](https://arxiv.org/abs/2408.13494v2)
