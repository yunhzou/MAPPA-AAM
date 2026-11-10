# Evidence and claim boundaries

Current decoded reaction recovery is 1,833/1,851 (one-seed Golden), 1,842/1,851
(ten-seed Golden) and 47,796/47,974 (one-seed EnzymeMap). All three reaction-level evaluations are resolved. Catalogue recovery is separately 1,834, 1,843 and
47,798. Golden catalogue verification permits endpoint chemical symmetry while
its decoded audit requires exact indices. EnzymeMap uses exact indices for both.
All required mapping and decoding outputs are available.

Successful compute CPU means for matching plus reduction / decoding are
1.31 / 1.88, 11.64 / 9.28 and 4.41 / 8.75 seconds per reaction, respectively.
Both directions are included. These records combine the documented matching and
postprocessing passes and successful recoveries; they are not an uninterrupted
end-to-end run of one final performance revision. Publication, failed attempts
and external verification are separately accounted for in the source reports.

The separately paired current Golden study records MAPPA catalogue recovery of
1,839/1,851, expanded-SLAP recovery of 1,794/1,851 and standard binary bidirectional
SLAP recovery of 1,620/1,851. The shared 1,834-reaction completed cohort has search
means of 1.202, 12.475 and 0.111 CPU seconds, respectively; MAPPA decoding adds
1.617 seconds. These support a speed advantage over expanded SLAP, not default
SLAP. The broader default-SLAP/MAPPA cohort and full-denominator outcomes remain
separately reported. Do not pool this study with the saved seed comparison.

At the saved one-seed detected minimum, 1,829 reactions have certified distinct
heavy-mapping counts; 22 counts remain incomplete. The certified-subset mean is
1.36, while the full-dataset median is one and its 95th percentile is three.
Minimum-set reference inclusion is 1,691 confirmed, 158 absent and two unresolved.
The expected accuracy of uniform selection is bounded at 78.31--79.45%; this is
not a measured chemistry-aware top-1 predictor. Hydrogen choices do not multiply
the counts, but the ranking score remains the saved all-atom model count.

The event-offset snapshot accounts for every processed reaction: 49,509 exact
offsets, 120 bounded intervals and 196 absent references. Figure 7 places each
bounded case at its saved lower endpoint, following the author's requested
summary. These values remain lower bounds in the source evidence. Retained-branch
offsets use the same all-atom decoder model on both sides; GT selects a branch
rather than supplying its count.

The same bars include anchored completion of all 196 absent
references: 18 Golden and 178 EnzymeMap. Their offset is the attained
reference-conditioned count minus the original detected all-atom minimum, not
the minimum of an original retained branch. All 392 directions preserve supplied
pairs and pass direct physical event counting. Nine EnzymeMap completions lie
two events below the original detected minimum. Eleven EnzymeMap references
require endpoint-index translation; partial Golden completions and index
translations are recorded in the report. These results do not increase the
independent recovery figures above or resolve the 120 bounded branch offsets.
The two merged distributions count all 1,851 Golden and 47,974 EnzymeMap reactions
exactly once, with one color per dataset and no separate completion legend.
Negative ticks are omitted when the dataset has no negative offsets.

The model retains a fixed ordinary-boundary approximation and selected starting
orbits. Positive decoded witnesses are replayed and counted, but model-window
coverage is not physical event-bin accuracy. EnzymeMap records 2,801 positive
direction witnesses whose model and physical counts differ.

Archived September comparator and coordinate evidence remains unchanged and is
labeled in the supplement. The 140-coordinate collection is not an accuracy
benchmark. Its 120,052 families and 300 event classes belong to that archived
study and must not be substituted for current decoded results.
