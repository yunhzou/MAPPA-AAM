# Published pathways for the gold rearrangement

Original pathway schemes from González Pérez et al., *Mechanism of the Gold-Catalyzed Rearrangement of (3-Acyloxyprop-1-ynyl)oxiranes: A Dual Role of the Catalyst*, **J. Org. Chem. 2009, 74, 2982–2991**. [Original paper · DOI: 10.1021/jo802516k](https://pubs.acs.org/doi/10.1021/jo802516k).

Compare these schemes with the [3D growth animation and its two oxygen assignments](../README.md). **O*** in the animation denotes the original epoxide oxygen; it is our tracking label, not an isotope label in the paper.

## Route a: initial 1,3-ester migration

![Original Scheme 2: route a through initial 1,3-ester migration](scheme-2-route-a.png)

**Animation Candidate 1:** the original epoxide oxygen becomes the product's ester-link oxygen. This endpoint correspondence is consistent with route a, which the paper favors energetically.

## Route b: initial 1,2-ester migration

![Original Scheme 4: route b through initial 1,2-ester migration](scheme-4-route-b.png)

**Animation Candidate 2:** the original epoxide oxygen becomes the product's ketone oxygen. Route b proceeds through a carbenoid and reaches intermediate 14.

## Route c: initial oxirane activation

![Original Scheme 5: route c through initial oxirane activation](scheme-5-route-c.png)

Route c also matches **Candidate 2**. It proceeds through an allene and converges with route b at intermediate 14. The two routes therefore share an endpoint oxygen pattern; endpoint AAM alone cannot distinguish them.

## Reading the comparison

The original schemes use **L = PH₃** and include literature energies. The animation uses the supplied **AuPPh₃ complexes 2 → 9**, with all 65 atoms included in the search. The energies shown above are from the paper, not MAPPA calculations. MAPPA recovers endpoint correspondences compatible with these hypotheses; it does not establish the intermediate steps or their kinetic feasibility.

The images are direct crops of Schemes 2, 4 and 5 on journal pages 2985–2987. Original structures, captions and footnotes are retained. These third-party figures remain attributable to the original publication and are not covered by this repository's code license. The complete PDF is not redistributed.

[Extraction provenance](provenance.json) records source and image hashes, PDF pages and crop regions. With the source PDF and Poppler installed, regenerate them from the repository root:

```bash
python manuscript/scripts/gold_film/extract_schemes.py \
  supplied-paper.pdf manuscript/animations/gold_rearrangement/pathway-reference
```
