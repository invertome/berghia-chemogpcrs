# Methods — Non-chemoreceptor GPCR classification

(For inclusion in the paper Methods, between the orthogroup-clustering description and the candidate-ranking description.)

## Three-source consensus classification of candidate GPCRs

To distinguish bona fide chemoreceptor candidates from GPCRs of other functional families (bioamine, neuropeptide, opsin, lipid, nucleotide, glycoprotein-hormone, class B/C/F), we developed a three-source consensus classifier that compares each candidate to a curated reference set of 728 reviewed GPCRs from UniProt/Swiss-Prot (human, mouse, rat, *Drosophila melanogaster*, *Caenorhabditis elegans*) spanning ten coarse families. We deliberately restricted references to manually-reviewed entries: NCBI XP_* and TrEMBL annotations propagate misidentification on invertebrate receptors and would have introduced systematic bias.

Each candidate is classified independently by three methods:

1. **HMM scan.** We trained a profile-HMM library (HMMER v3.4; Eddy 2011) from the curated reference alignments — 8 family-level HMMs and 11 subfamily-level HMMs for the two largest families (aminergic: 5HT, dopamine, histamine, norepinephrine, octopamine, tyramine; peptide: NPY-NPF, cholecystokinin, kinin, melatonin, opioid, somatostatin, tachykinin, vasopressin-oxytocin) plus four Pfam fallback profiles (PF00001 7tm_1, PF00002 7tm_2, PF00003 7tm_3, PF01534 Frizzled). Per-family E-value thresholds were calibrated by leave-one-out cross-validation against both ≥0.90 recall and ≥0.95 precision targets; all retained families exceeded both targets, with five at perfect 1.000/1.000.

2. **Orthogroup vote.** Candidates inheriting an OrthoFinder orthogroup whose members carry concordant Swiss-Prot or HMM-derived family annotations (≥80 % consensus among ≥3 annotated members, candidate excluded from its own vote) are classified into that family.

3. **Phylogenetic placement.** Candidates are placed onto reference trees using EPA-ng (Barbera et al. 2019). The pipeline maintains a hierarchical tree set: a coarse backbone tree (54 leaves, 6 representatives × 9 families) plus medium-granularity subtrees for the two largest non-chemoreceptor families (aminergic, 165 leaves; peptide, 189 leaves). Reference-tree alignments use the same filter stack that we apply to the main pipeline (Whelan et al. 2018 — PREQUAL residue pre-mask; Wheeler 2025 — CLOAK consensus mask over a five-MAFFT-variant ensemble; Zhang et al. 2021 — TAPER residue-outlier mask; ClipKit column trim; Steenwyk et al. 2020). Trees were inferred with IQ-TREE 3 (Nguyen et al. 2015; Minh et al. 2020) under ModelFinder Plus (Kalyaanamoorthy et al. 2017), with 1000 ultrafast bootstrap replicates, 1000 SH-aLRT replicates, and transfer bootstrap expectation (Lemoine et al. 2018). Placement calls require a likelihood-weight ratio (LWR) ≥ 0.80 on the best edge.

A consensus call is produced from these three sources: 3-of-3 agreement is recorded as **high-confidence non-chemoreceptor**; 2-of-3 agreement is recorded as **medium-confidence (likely non-chemoreceptor)**; otherwise the candidate is retained as a **chemoreceptor candidate**. This stratification is conservative by design — false-positive non-chemoreceptor calls (genuine chemoreceptors mis-flagged) would prune real signal from the wet-lab triage list, whereas false-negative non-chemoreceptor calls (a 5HT receptor surviving as a "candidate") are recoverable downstream by the orthogonal evidence channels (synteny, dN/dS, expression, paralog tandem-cluster context).

We retrospectively validated the classifier against a curated *Aplysia californica* truth set drawn from Cummins et al. (2009, *J. Comp. Physiol. A* 195:867 — three published rhinophore chemoreceptors NP_001191513.1/NP_001191494.1/NP_001191495.1) plus four published non-chemoreceptors (5HT-1, 5HT-2, 5HT4, allatostatin-A). Reported metrics include chemoreceptor correctness (true-negative rate among known chemoreceptors), non-chemoreceptor recall, family accuracy, and subfamily accuracy.

Implementation: `scripts/classify_via_hmm.py`, `scripts/classify_via_og_vote.py`, `scripts/classify_via_placement.py`, `scripts/classify_consensus.py`, orchestrated by `06c_classify_non_chemoreceptors.sh`. Reference HMMs and trees in `results/classification/`.

## Caveats and design choices

The medium-granularity subfamily-level resolution is currently restricted to aminergic and peptide families because (a) these account for ~50 % of the curated reference set and are the most common non-chemoreceptor noise in invertebrate chemoreceptor candidate pools, and (b) other coarse families (class-B-secretin, lipid, class-C, class-F-frizzled, glycoprotein-hormone, nucleotide, opsin) lacked sufficient subfamily-balanced reference depth (n ≥ 10 per subfamily) for stable HMM training. Coarse-family classification is available for all ten families.

Family monophyly on the backbone tree is recovered for 6/9 coarse families (aminergic, class-B-secretin, class-C, class-F-frizzled, glycoprotein-hormone, opsin); lipid, nucleotide, and peptide are non-monophyletic, consistent with their published phylogenies as functional groupings spanning multiple evolutionary origins (e.g., cannabinoid, leukotriene, and prostaglandin receptors arise from independent class-A lineages). Within the aminergic and peptide subtrees, subfamily-level monophyly is incomplete — a known feature of receptor evolution where adrenergic, dopaminergic, and trace-amine clades remain interleaved. EPA-ng placement uses local edge likelihoods, so global non-monophyly does not invalidate placement; the 3-source consensus design absorbs residual disagreement by capping confidence to "medium" when sources disagree.

## References to add to the bibliography

- Eddy SR. 2011. Accelerated profile HMM searches. *PLoS Comput Biol* 7(10):e1002195.
- Barbera P, Kozlov AM, Czech L, Morel B, Darriba D, Flouri T, Stamatakis A. 2019. EPA-ng: massively parallel evolutionary placement of genetic sequences. *Syst Biol* 68(2):365–369.
- Whelan S, Irisarri I, Burki F. 2018. PREQUAL: detecting non-homologous characters in sets of unaligned homologous sequences. *Bioinformatics* 34(22):3929–3930.
- Wheeler TJ. 2025. CLOAK: consensus-based masking for multiple sequence alignment uncertainty. *bioRxiv* 691663.
- Zhang C, Zhao Y, Braun EL, Mirarab S. 2021. TAPER: pinpointing errors in multiple sequence alignments despite varying rates of evolution. *Methods Ecol Evol* 13(1):91–105.
- Steenwyk JL, Buida III TJ, Li Y, Shen X-X, Rokas A. 2020. ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. *PLoS Biol* 18(12):e3001007.
- Nguyen L-T, Schmidt HA, von Haeseler A, Minh BQ. 2015. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Mol Biol Evol* 32(1):268–274.
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol* 37(5):1530–1534.
- Kalyaanamoorthy S, Chernomor O, von Haeseler A, Minh BQ. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods* 14(6):587–589.
- Lemoine F, Domelevo Entfellner J-B, Wilkinson E, Correia D, Dávila Felipe M, De Oliveira T, Gascuel O. 2018. Renewing Felsenstein's phylogenetic bootstrap in the era of big data. *Nature* 556:452–456.
- Cummins SF, Erpenbeck D, Zou Z, Claudianos C, Moroz LL, Nef P, Degnan BM. 2009. Candidate chemoreceptor subfamilies differentially expressed in the chemosensory organs of the mollusc Aplysia. *BMC Biol* 7:28.
