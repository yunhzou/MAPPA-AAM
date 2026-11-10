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

At the saved one-seed detected minimum, all 1,851 Golden reactions have recorded
qualifying branch compositions. Median 31, mean 51.06, 95th percentile 120,
maximum 2,250. Internal shuffles and repeated histories count once. This is a
count of fragment atom-set compositions, not chemically distinct full mappings.
The analysis makes no new symmetry queries or optimization calls. Figure 8 uses
`minimum_branches.json` and `minimum_branches_reactions.json.gz`; older
`minimum_choices` files remain historical evidence and no longer supply it.

Figure 7 uses `reference_offsets.json.gz`: a verified GT-consistent mapping's
event count minus the saved detected minimum, under the same all-atom decoder
model. Select the lowest recorded verified count across both directions without
new search. All 49,825 offsets are differences of concrete saved counts. The
superseded `branch_offsets.json.gz` estimated minima of reference-containing
branches and is not the plotted metric.

Saved anchored completions supply the 18 Golden and 178 EnzymeMap references
absent from the original decoded output. All 392 directions preserve supplied
pairs and pass direct physical event counting. Nine EnzymeMap completions lie
two events below the original detected minimum. Eleven EnzymeMap references
require endpoint-index translation; partial Golden completions and translations
remain documented in the report. These results do not increase independent
recovery. The two distributions count every reaction once, with one color per
dataset, no separate completion legend and no unoccupied negative ticks.

The corrected figure moves 126 Golden and 2,196 EnzymeMap reactions to higher
offsets than the previous branch-minimum figure. Counts at zero are 1,630/1,851
(88.06%) and 44,086/47,974 (91.90%). A saved hydrogen completion is not necessarily
optimal, so this measurement is distinct from reference membership among all
minimum-event choices.

The model retains a fixed ordinary-boundary approximation and selected starting
orbits. Positive decoded witnesses are replayed and counted, but model-window
coverage is not physical event-bin accuracy. EnzymeMap records 2,801 positive
direction witnesses whose model and physical counts differ.

Archived September comparator and coordinate evidence remains unchanged and is
labeled in the supplement. The 140-coordinate collection is not an accuracy
benchmark. Its 120,052 families and 300 event classes belong to that archived
study and must not be substituted for current decoded results.
