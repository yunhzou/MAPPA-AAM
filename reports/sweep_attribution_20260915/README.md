# SLAP literature score and cut-sweep attribution

Original SLAP source: Koda and Saito, DOI https://doi.org/10.26434/chemrxiv-2025-hthwn . The ChemRxiv manuscript, p. 12, identifies the 1,758-reaction LocalMapper Golden version. Figure 4 (p. 15) reports 86.9% reference inclusion among returned SLAP solutions. Verified against the primary PDF indexed text on 2026-09-15; direct PDF retrieval returned 403. No inference of a same-protocol gain from the literature result is made.

Our strict 1,851-record default binary SLAP result is 1,591 (85.95%). Paired uncut binary/weighted bidirectional SLAP recovers 1,661 (89.74%); our sweep raises it to 1,796 (97.03%), adding 135 references. MAPPA one-seed uncut recovers 1,489 (80.44%); default sweep recovers 1,834 (99.08%), adding 345. Existing evidence/unswept.json and competitors.json remain unchanged.

Introduction and Methods now credit the cut-sweep strategy as introduced in this work. The baseline protocol, results, and comparison captions distinguish original SLAP from SLAP augmented with our sweep. No algorithm changes or experiment reruns.
