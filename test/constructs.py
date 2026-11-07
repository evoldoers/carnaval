"""Reference constructs and their target secondary structures.

Each construct names the helices whose formation the folding-time measurement
waits for.  Pair indices are zero-based positions in the sequence.
"""

# ---------------------------------------------------------------------------
# Yeast phenylalanine tRNA, PDB entry 1ehz.
#
# The sequence and the reference secondary structure are the ones used in the
# CARNAVAL paper (carnaval-paper/data/1ehz.seq, whose second line is the
# structure annotated from the PDB entry).  Parsing that annotation gives four
# canonical helices plus one tertiary base pair; the tertiary pair is recorded
# but is not required for the "folded" call, since a single pair frays too
# readily to serve as a landmark.
# ---------------------------------------------------------------------------
TRNA = {
    "name": "trna-1ehz",
    "seq": "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
    "board": (128, 128, 1),
    "helices": {
        "acceptor":  [(0, 71), (1, 70), (2, 69), (3, 68), (4, 67), (5, 66), (6, 65)],
        "D-arm":     [(9, 24), (10, 23), (11, 22), (12, 21)],
        "anticodon": [(26, 42), (27, 41), (28, 40), (29, 39)],
        "T-arm":     [(48, 64), (49, 63), (50, 62), (51, 61), (52, 60)],
    },
    "tertiary": [(17, 55)],
    "core": [],
}

# ---------------------------------------------------------------------------
# Hammerhead ribozyme.
#
# This is a hammerhead-topology construct, not a natural isolate: a three-way
# junction carrying the two published core consensus motifs as its junction
# linkers (5'-CUGANGA with N=U, and 5'-GAAA), a GNRA tetraloop capping stem II
# as in the viral hammerheads, and a GUC cleavage triplet in the stem III loop.
# The three stems are designed, so the target structure is known by
# construction rather than taken from an annotation.
#
# Substituting a natural sequence is a matter of replacing "seq" and the helix
# pair lists, and wants an annotated structure for that isolate to do properly.
# ---------------------------------------------------------------------------
_HH_SEGMENTS = [
    ("stemI-5",   "GGAGGC"),   #  0-5
    ("core-1",    "cugauga"),  #  6-12   CUGANGA consensus, N=U
    ("stemII-5",  "GGAC"),     # 13-16
    ("loopII",    "guga"),     # 17-20   GNRA tetraloop
    ("stemII-3",  "GUCC"),     # 21-24
    ("core-2",    "gaaa"),     # 25-28   GAAA consensus
    ("stemIII-5", "AGGCAG"),   # 29-34
    ("loopIII",   "guca"),     # 35-38   carries the GUC cleavage triplet
    ("stemIII-3", "CUGCCU"),   # 39-44
    ("stemI-3",   "GCCUCC"),   # 45-50
]

HAMMERHEAD = {
    "name": "hammerhead",
    "seq": "".join(s for _, s in _HH_SEGMENTS),
    "board": (64, 64, 1),
    "helices": {
        "stemI":   [(k, 50 - k) for k in range(6)],
        "stemII":  [(13 + k, 24 - k) for k in range(4)],
        "stemIII": [(29 + k, 44 - k) for k in range(6)],
    },
    "tertiary": [],
    # The core must be unpaired for the junction to be catalytically competent.
    "core": list(range(6, 13)) + list(range(25, 29)),
}

CONSTRUCTS = {c["name"]: c for c in (TRNA, HAMMERHEAD)}


def check():
    """Verify every declared pair is complementary and well formed."""
    comp = {frozenset("gc"), frozenset("au"), frozenset("gu")}
    for c in CONSTRUCTS.values():
        seq = c["seq"].lower()
        allpairs = [p for h in c["helices"].values() for p in h] + c["tertiary"]
        for i, j in allpairs:
            assert 0 <= i < j < len(seq), (c["name"], i, j, len(seq))
            assert frozenset((seq[i], seq[j])) in comp, (c["name"], i, j, seq[i], seq[j])
            # CARNAVAL forbids pairing between units closer than four along the chain
            assert j - i >= 4, (c["name"], i, j)
        seen = set()
        for i, j in allpairs:
            assert i not in seen and j not in seen, (c["name"], i, j)
            seen.add(i); seen.add(j)
    return True


if __name__ == "__main__":
    check()
    for c in CONSTRUCTS.values():
        npair = sum(len(h) for h in c["helices"].values())
        print(f"{c['name']:12s} {len(c['seq']):3d} nt  "
              f"{len(c['helices'])} helices / {npair} bp  "
              f"helices: {', '.join(c['helices'])}")
        print(f"  {c['seq']}")
