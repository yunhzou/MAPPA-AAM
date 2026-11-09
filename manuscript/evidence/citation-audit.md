# Initial bibliography audit

Checked against primary publication records or author-maintained resources
on 10 September 2026. These are the five references used in this first draft.

| BibTeX key | Claim supported | Primary record |
| --- | --- | --- |
| `lin2022` | Golden benchmark provenance and common evaluation concerns | [PubMed](https://pubmed.ncbi.nlm.nih.gov/34726834/), [authors' data repository](https://github.com/Laboratoire-de-Chemoinformatique/Reaction_Data_Cleaning); DOI 10.1002/minf.202100138 |
| `schwaller2021` | Unsupervised reaction-language-model attention for mapping | [Full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC8026122/); DOI 10.1126/sciadv.abe4166 |
| `chen2024` | Human-in-the-loop LocalMapper | [Nature Communications](https://www.nature.com/articles/s41467-024-46364-y); DOI 10.1038/s41467-024-46364-y |
| `koda2025` | SLAP's WL-like labels, sequential assignment and alternatives | [ChemRxiv version 1](https://www.cambridge.org/engage/chemrxiv/article-details/691d9f5ea10c9f5ca1852be4), [authors' implementation](https://github.com/shin1koda/slap-mapper); DOI 10.26434/chemrxiv-2025-hthwn |
| `mckay2014` | Canonical labeling / nauty | [Author's publication list](https://users.cecs.anu.edu.au/~bdm/publications.html), [author preprint](https://arxiv.org/abs/1301.1493); DOI 10.1016/j.jsc.2013.09.003 |

Updated 8 November 2026: the primary SLAP preprint text confirms that its Golden
evaluation uses the 1,758-reaction LocalMapper version and counts reference
inclusion among returned alternatives. Its reported 86.9% is not first-returned
accuracy. The paper now distinguishes this literature result from our
1,851-reaction evaluation, our measured first-returned score, and our expanded
SLAP wrapper. The tested software is pinned separately to 1.0.0 / `ea248fd`.
The independent returned-output audit and protocol are in `paired_slap.json`;
the expanded sweep is our experiment and is not attributed to the SLAP authors.

RXNMapper and LocalMapper support the introduction's selected-answer objective;
SLAP supports the alternative-generating objective. This is a taxonomy of output
objectives, not a claim that learned and optimization-based methods form disjoint
algorithm classes. The proposed downstream chemical selector is not evaluated.
