# Fresh paired Golden comparison

Status: preparing the pilot. No current-version speedup is claimed yet.

The comparison runs all 1,851 Golden reactions with both directions on one CPU
model in one allocation. Current MAPPA uses tolerance 2, one seed, cap 100 and
the complete uncut plus single-edge sweep. Its indexed fragment cache is enabled
and starts empty for each direction, sharing discoveries across that direction's
cuts. Decoding uses the original weights, a five-event window and the existing
five-second unique-fragment allowance.

SLAP uses upstream revision `ea248fd9494f52f4865193e87a98cc92c62b5f9e`, both
binary and weighted modes, heavy-atom symmetry breaking, explicit hydrogens and
the same cut scope. Its isolated element-balancing placeholders preserve the
upstream input policy. Original prepared graph indices and weights are shared
between the methods; reference labels never enter matching.

Search time includes graph preparation, matching and reduction. Input loading,
publication, MAPPA decoding and external verification are separate. Numerical
libraries use one thread each. Common completed reaction cohorts, tail times,
resource failures and reference coverage will accompany the mean CPU comparison.
Catalogue reference inclusion uses endpoint chemical equivalence for both
methods. MAPPA's exact-index decoded-window audit is a separate metric.

The maintained runner is `bench/compare_slap.py`, using the existing persistent
pool, aggregate storage and ledger. It stops admission on an execution/resource
failure and drains active work before investigation. Completed outputs and failed
attempts remain saved; there are no automatic retries.

Six local tests passed: upstream graph/output parity for balanced, aromatic and
unbalanced inputs; rejection of invalid repeated heavy labels; joint connectivity
and unmatched-atom reference checks; and both methods' output publication followed
by MAPPA decoding and external verification. No core algorithm was modified.
