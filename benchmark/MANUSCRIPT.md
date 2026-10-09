**Host-resolved reverse vaccinology: VaxForge couples whole-proteome antigen discovery to species-specific epitope prediction and construct assembly**

Ali Arslan

Department of Molecular Biology and Genetics, Institute of Natural Sciences, Yıldız Technical University, Istanbul 34220, Turkey

Correspondence: aliarslan47@gmail.com

Draft restructured 2026-10-10 (supersedes the 2026-09-20 draft). Extended methods narrative in `MANUSCRIPT_METHODS.md`; ground-truth antigen evidence and DOIs in `PROTECTIVE_ANTIGENS.md`; accessions in `PATHOGEN_ACCESSIONS.md`; per-antigen miss attribution in `MISS_ATTRIBUTION.md`. Figure file mapping: Fig. 1 `fig1_antigen_recovery`, Fig. 2 `fig2_recall_bar`, Fig. 3 `fig3_fold_vs_recall`, Fig. 4 `fig4_bacteria_panel`, Fig. 5 `fig6_virus_antigen_rank`, Fig. 6 `fig7_virus_multihost`, Fig. S1 `figS1_brucella_vaccine`, Fig. S2 `fig5_multihost_jaccard`.

---

**Abstract**

Reverse vaccinology ranks candidate antigens from a pathogen genome, but the tools in common use stop at the ranking step and model the host as a single human reference, which limits their value for the veterinary and zoonotic targets that carry much of the global vaccine burden. We present VaxForge, a deterministic pipeline that carries a proteome through a soft-filter antigen funnel, class I and class II epitope prediction against curated allele panels for four mammalian hosts, multi-epitope construct assembly and an in silico immunisation readout. Benchmarked under default parameters against NERVE 2.0 and Vaxign-ML on ten bacteria at whole-proteome scale, VaxForge recovered 37 of 48 independently validated protective antigens, against 47 of 48 for Vaxign-ML and 14 of 48 for NERVE 2.0, and recovered all four licensed 4CMenB antigens where NERVE 2.0 recovered one. We then attributed every missed antigen to a specific pipeline step. Six of the eleven misses were removed by the mandatory host-homology filter, a relationship that is absolute in one direction: no recovered antigen carried a human homologue (0 of 37; Fisher exact p = 3.8 × 10⁻⁵). The recall gap is therefore a quantified safety-sensitivity trade-off rather than a failure of discrimination. Across ten viruses the predicted epitope repertoire of the protective antigen differed between host species in every case, with host pairs sharing no epitopes at all for several viruses, a dimension neither comparator can resolve.

---

**Introduction**

Reverse vaccinology began with serogroup B *Neisseria meningitidis*, where candidate antigens were derived from the whole-genome sequence rather than from cultured components (Pizza et al., 2000), and produced the antigens behind the licensed four-component meningococcal vaccine: factor H binding protein, Neisserial adhesin A and Neisserial Heparin Binding Antigen, together with a PorA-bearing outer-membrane vesicle (Giuliani et al., 2006). Two decades later the approach is served by mature software spanning two design philosophies (Goodswen et al., 2023). Filter-based environments such as NERVE 2.0 pass proteins through a cascade of cut-offs on localisation, topology and adhesin probability (Conte et al., 2024). Classifier-based methods score sequence directly: Vaxign-ML applies gradient boosting to curated descriptors (Ong et al., 2020), and the current generation replaces those descriptors with protein language model embeddings, as in Vaxign-DL (Chen et al., 2024), the descriptor-free ESM-2 pipeline of Podda et al. (2025), and the sequence-plus-structure framework PLGDL, which added geometric deep learning and validated a previously uncharacterised orthopoxvirus antigen in vivo (Zai et al., 2025).

That progress has been concentrated almost entirely on one question: how accurately a model separates protective antigens from the rest of a proteome. Three practical gaps have been left largely untouched. First, a hard virulence or localisation filter can raise the apparent enrichment of a candidate list while silently discarding true antigens, and the Neisserial Heparin Binding Antigen is the standing example, being absent from curated virulence databases yet a licensed vaccine component. Second, epitope-level output is computed against a single human reference, which leaves cattle, pig and mouse targets without an honest presentation model even though NetMHCpan and NetMHCIIpan support many non-human alleles (Reynisson et al., 2020). Third, the pipelines end at a ranked protein list, so epitope selection, construct assembly and a first immunological check remain manual. A classifier that ranks antigens more accurately does not close any of the three.

VaxForge was built to close them inside one deterministic pipeline. It replaces the hard virulence gate with a soft, score-weighted term, so a genuine surface antigen with no virulence annotation can still pass; it predicts class I and class II binding against curated allele panels for human, mouse, bovine and porcine hosts; and it continues past ranking into a characterised multi-epitope construct and an in silico immunisation readout. Here we quantify its behaviour against NERVE 2.0 and Vaxign-ML on ten bacteria and ten viruses under default parameters, we attribute every antigen it missed to a specific pipeline step rather than reporting recall alone, and we show that the resulting recall gap against a pure classifier is dominated by a deliberate safety filter whose action we can measure.

---

**Results**

**VaxForge couples a soft-filter funnel to a host-resolved immunological layer.** VaxForge detects whether its input is a proteome, coding sequences or a raw genome, then runs discovery, an antigen funnel, host-aware epitope prediction, construct assembly and an immunisation check, reading every cutoff from one configuration file so that no threshold is hard-coded. The funnel applies three hard filters, on subcellular localisation, transmembrane topology and host homology, and treats virulence evidence, signal peptides, antigenicity, cross-strain conservation and adhesin probability as soft, score-weighted terms. Surviving proteins enter the immunological layer, which is the part that distinguishes the pipeline: class I and class II binding are predicted for the alleles of the chosen host, so a T-cell epitope is retained only if it is presented in the intended species, and three further predictors refine the readout through class I immunogenicity, proteasomal cleavage with TAP transport, and linear B-cell epitopes. Retained epitopes are assembled with defined linkers and an adjuvant into a construct that is characterised physicochemically and passed to the immunisation step. Population coverage is computed for human HLA and reported as unavailable, rather than imputed, for hosts without allele-frequency data. Table 1 sets the resulting capability surface beside the two comparators.

**Table 1.** Inputs, outputs and capabilities of the three compared tools (✓ present, ✗ absent). Comparator capabilities are stated as described in their source publications (NERVE 2.0, Conte et al., 2024; Vaxign-ML, Ong et al., 2020).

| Capability | VaxForge | NERVE 2.0 | Vaxign-ML |
|---|---|---|---|
| Input | proteome, CDS or genome (auto-detected) | bacterial proteome + Gram type | bacterial protein proteome |
| Supported pathogens | bacteria and viruses | bacteria only | bacteria (viruses not described in publication) |
| Approach | rule + ML hybrid funnel with an immunological layer | rule-based modular pipeline with AI modules | supervised ML classifier (XGBoost) |
| Antigen ranking | ✓ weighted candidacy score | ✓ Select + score | ✓ protegenicity score |
| Autoimmunity (host-homology) filter | ✓ human, mandatory | ✓ human + optional mouse | ✗ |
| Retains virulence-negative antigens | ✓ soft funnel | ✗ hard Select removes them | ✓ score-based |
| T-cell / MHC epitope prediction | ✓ NetMHCpan/II + processing + immunogenicity | ✓ MHC-I/II, linear, human HLA supertypes | ✗ |
| B-cell epitope prediction | ✓ BepiPred-3.0 | ✗ | ✗ |
| Host-specific MHC allele selection | ✓ human/mouse/cattle/pig | partial (epitopes human HLA only; mouse for safety) | ✗ |
| Population coverage | ✓ human HLA (other hosts reported "no data") | ✓ human HLA | ✗ |
| Multi-epitope construct | ✓ linkers + adjuvant | ✗ | ✗ |
| In silico immunisation | ✓ | ✗ | ✗ |
| Primary output | ranked antigens + host-specific epitopes + construct + immunisation | PVC shortlist + T-cell epitopes / population coverage | protegenicity score per protein |

**VaxForge recovers all four licensed meningococcal antigens where NERVE 2.0 recovers one.** We first benchmarked on *N. meningitidis* MC58 with the four 4CMenB antigens as ground truth. On an expanded 204-protein set (four protective plus 200 background, fixed seed) VaxForge recovered all four (68 candidates; fold enrichment 3.12; Fisher exact p = 0.012), ranking PorA first, factor H binding protein third, the Heparin Binding Antigen sixth and adhesin A eighth. Vaxign-ML also recovered all four from 59 candidates, whereas NERVE 2.0 recovered only PorA. The pattern held on a compact 64-protein set and on the full 2,003-protein proteome for the two faster tools (Table 2; Figs. 1 and 2). The decisive case is the Heparin Binding Antigen, which is virulence-database negative: VaxForge retained it through the soft funnel while NERVE 2.0 removed it at its hard selection step. That an independent classifier reached the same 4 of 4 indicates a property of soft filtering rather than a single-tool artefact.

**Table 2.** Meningococcal head-to-head across three input scales (default parameters, identical input).

| Set | Tool | Recall | Candidates | Fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| 204 | VaxForge | 4/4 | 68 | 3.12 | 0.012 | 0.99 |
| 204 | Vaxign-ML | 4/4 | 59 | 3.64 | 0.007 | 0.97 |
| 204 | NERVE 2.0 | 1/4 | 3 | 25.0 | 0.058 (ns) | n/a |
| Full 2,003 | Vaxign-ML | 4/4 | 613 | 3.28 | n/a | n/a |
| Full 2,003 | NERVE 2.0 | 1/4 | 9 | 62.5 | n/a | n/a |

![Figure 1](figures/fig1_antigen_recovery.png)

**Figure 1.** VaxForge and Vaxign-ML recover all four licensed meningococcal antigens, whereas NERVE 2.0 recovers only PorA. Each point is one tool's result for one antigen, and colour identifies the tool (legend, top right). The number beside each point is that antigen's rank in the tool's candidate list, meaning its position from the top (1 is first, and a lower number is better). A cross in the "not recovered" column marks an antigen the tool did not return at all. Ranks and recall were computed under default parameters on identical input.

![Figure 2](figures/fig2_recall_bar.png)

**Figure 2.** Recall over known antigens separates the tools more clearly than area-under-curve. Per-tool recall of the four meningococcal antigens with secondary AUC annotation; AUC carries a small-sample caveat and is not used as the headline metric (n = 4 positives).

**Fold enrichment and recall diverge, and can move in opposite directions.** NERVE 2.0 produced short candidate lists with high fold enrichment, 25 on the 204 set and 62.5 on the full proteome, consistent with the value near 16.6 reported for its best configuration (Conte et al., 2024), yet its Fisher exact p was not significant because the candidate count was very small, and it missed three of the four licensed antigens at every scale (Fig. 3). Those three losses were not caused by low continuous scores from NERVE 2.0, which rated them highly, but by the hard selection step applied after scoring: periplasmic localisation, an adhesin probability below the 0.5 cut-off, and an unknown-localisation call with one predicted transmembrane helix. Fold enrichment rewards selectivity and can therefore rise precisely as recall falls. The quantity that reflects what a tool delivers to a user is recall over independently validated antigens, together with the ranks at which they appear.

![Figure 3](figures/fig3_fold_vs_recall.png)

**Figure 3.** High fold enrichment does not imply high recall. Each bar is the fold enrichment of the known protective antigens in one tool's candidate list, that is, how many times more concentrated the known antigens are among the candidates than expected by chance. The label above each bar is that run's recall, the count of known antigens the tool actually recovered out of the total (4/4 means all four were recovered, 1/4 only one). The dashed line marks the published NERVE 2.0 fold value of 16.6. NERVE 2.0 reaches the highest fold with the lowest recall.

**The pattern holds across ten bacteria at whole-proteome scale.** To move beyond a single exemplar we ran nine further zoonotic bacteria at full proteome scale, each a natural pathogen of all four mammalian hosts, chosen to span genome size from 2.0 to 5.5 Mb and to balance Gram type and cell biology, from the obligate intracellular *Coxiella burnetii*, through the spirochaete *Leptospira interrogans*, to the plasmid-encoded toxins of *Bacillus anthracis* (Table 3). Across the 44 protective antigens of these nine species VaxForge recovered 33, against 43 for Vaxign-ML and 13 for NERVE 2.0; including the meningococcal exemplar the totals are 37, 47 and 14 of 48 (Fig. 4). VaxForge reached 4 of 5 or 5 of 5 on seven of the nine species, including *Staphylococcus aureus* (iron-regulated surface determinant B at rank 4 of 2,889 proteins), *Salmonella* Typhimurium (5 of 5), *Brucella suis*, *Mycobacterium bovis*, *L. interrogans* and *B. anthracis*, for which the anthrax protective antigen was recovered at rank 2. NERVE 2.0 returned none of the protective antigens for four of the nine species, including the two, *B. suis* and *B. anthracis*, where it reached the highest fold enrichment of any tool (185 and 274).

**Table 3.** Nine-bacterium head-to-head at full proteome scale (default parameters, identical input), ordered by genome size. Each species carries four or five protective antigens.

| Bacterium (genome, proteins) | Tool | Recall | Candidates | Fold | AUC |
|---|---|---|---|---|---|
| *C. burnetii* (2.0 Mb, 1,812) | VaxForge | 2/4 | 524 | 1.7 | n/a |
|  | Vaxign-ML | 4/4 | 477 | 3.8 | 0.91 |
|  | NERVE 2.0 | 0/4 | 9 | n/a | n/a |
| *S. suis* (2.2 Mb, 1,922) | VaxForge | 1/5 | 288 | 1.3 | n/a |
|  | Vaxign-ML | 5/5 | 628 | 3.1 | 0.88 |
|  | NERVE 2.0 | 0/5 | 0 | n/a | n/a |
| *S. aureus* (2.8 Mb, 2,889) | VaxForge | 5/5 | 729 | 3.98 | 1.00 |
|  | Vaxign-ML | 5/5 | 838 | 3.46 | 0.94 |
|  | NERVE 2.0 | 3/5 | 34 | 55.8 | n/a |
| *L. monocytogenes* (2.9 Mb, 2,844) | VaxForge | 4/5 | 478 | 4.79 | n/a |
|  | Vaxign-ML | 5/5 | 917 | 3.11 | 0.95 |
|  | NERVE 2.0 | 3/5 | 73 | 24.3 | n/a |
| *B. suis* (3.3 Mb, 3,249) | VaxForge | 4/5 | 929 | 2.8 | n/a |
|  | Vaxign-ML | 5/5 | 839 | 3.9 | 0.88 |
|  | NERVE 2.0 | 2/5 | 9 | 185.4 | n/a |
| *M. bovis* (4.4 Mb, 4,088) | VaxForge | 4/5 | 535 | 6.2 | n/a |
|  | Vaxign-ML | 4/5 | 695 | 4.7 | 0.95 |
|  | NERVE 2.0 | 0/5 | 12 | n/a | n/a |
| *L. interrogans* (4.7 Mb, 3,676) | VaxForge | 4/5 | 1,105 | 2.7 | n/a |
|  | Vaxign-ML | 5/5 | 1,211 | 3.0 | 0.88 |
|  | NERVE 2.0 | 0/5 | 116 | n/a | n/a |
| *S.* Typhimurium (4.9 Mb, 4,533) | VaxForge | 5/5 | 1,368 | 3.32 | 0.99 |
|  | Vaxign-ML | 5/5 | 1,276 | 3.56 | 0.92 |
|  | NERVE 2.0 | 3/5 | 118 | 23.6 | n/a |
| *B. anthracis* (5.5 Mb, 5,493) | VaxForge | 4/5 | 1,163 | 3.8 | n/a |
|  | Vaxign-ML | 5/5 | 1,297 | 4.2 | 0.94 |
|  | NERVE 2.0 | 2/5 | 10 | 274.4 | n/a |
| **Total (nine species)** | **VaxForge 33/44 · Vaxign-ML 43/44 · NERVE 2.0 13/44** |  |  |  |  |

![Figure 4](figures/fig4_bacteria_panel.png)

**Figure 4.** The recall pattern across ten bacteria at full proteome scale, ordered by genome size. Recall over the known protective antigens for each tool, shown as a fraction with the recovered-of-total count above each bar because the antigen count differs between species; the meningococcal exemplar is included alongside the nine generalisation species and an asterisk marks the six zoonotic bacteria added in this work. VaxForge matches Vaxign-ML on the species whose protective antigens are secreted or surface-exposed, and the species where it falls behind are those whose ground-truth sets are dominated by conserved housekeeping proteins; NERVE 2.0 recovers none of the antigens for four species.

**Every missed antigen is attributable to a named step, and the largest single cause is the host-homology filter acting as designed.** Reporting recall alone leaves the reader unable to distinguish a design choice from a defect, so we traced each of the 11 antigens that VaxForge did not recover to the step that removed it, using the pipeline's own intermediate outputs together with DIAMOND alignment against the human proteome (Methods; `MISS_ATTRIBUTION.md`). The misses partition cleanly into four causes, and none is a scale-dependent engineering failure of the kind corrected earlier in this work.

Six of the 11, the 6-phosphogluconate dehydrogenase of *S. suis*, the GroEL/HtpB chaperone of *C. burnetii*, enolase, glyceraldehyde-3-phosphate dehydrogenase, the Mip peptidyl-prolyl isomerase and the DnaJ chaperone, were removed by the mandatory host-homology filter, which excludes proteins above 35% identity to the human proteome. The association is absolute in one direction: all six carried a human homologue at 37.0 to 52.3% identity, while none of the 37 recovered antigens carried any human hit at all, most producing no alignment above threshold (0 of 37; Fisher exact p = 3.8 × 10⁻⁵, odds ratio infinite). All six belong to the same biological class, conserved cytoplasmic housekeeping proteins that moonlight on the bacterial surface, for which a detectable human orthologue is near-unavoidable. These are therefore not antigens the pipeline failed to recognise but antigens it declined on safety grounds, and they account for more than half the recall gap to a classifier that applies no host-homology exclusion.

Two further misses, the CapB capsule biosynthesis protein of *B. anthracis* and ActA of *L. monocytogenes*, expose a genuine configuration inconsistency rather than a design choice. The permitted-localisation list for bacteria was written with Gram-negative architecture in mind and admits extracellular, outer-membrane, cell-wall and periplasmic calls. Gram-positive bacteria have no outer membrane, so their surface-exposed proteins receive the PSORTb `membrane` label, which that list does not admit; both antigens were removed despite a maximal localisation score of 10.0 and a single predicted transmembrane helix. One miss, LigA of *L. interrogans*, follows an upstream misassignment, PSORTb calling a characterised surface lipoprotein cytoplasmic with a score of 9.97, which is a limit of the localisation predictor rather than of the funnel. The remaining two, ESAT-6 of *M. bovis* and L7/L12 of *Brucella*, survived all three hard filters and were lost at the epitope stage, having yielded no candidate peptide. We tested and rejected the obvious explanation: the toxicity filter flagged only 60% of ESAT-6 9-mers, the lowest of the four proteins examined, while the CFP-10 control was flagged at 82% and still produced two peptides, and protein length does not separate them (100 against 95 residues). The loss therefore occurs at the MHC binding and scoring thresholds, which is a legitimate prediction-level miss rather than a systematic fault.

**Table 4.** Attribution of all 11 protective antigens that VaxForge did not recover across the ten-bacterium benchmark, by the pipeline step responsible. Per-antigen human identity values and recovery status are in `results/miss_attribution.json`.

| Cause | n | Antigens | Interpretation |
|---|---|---|---|
| Host-homology hard filter (≥ 35% human identity) | 6 | 6PGD (*S. suis*), GroEL/HtpB (*C. burnetii*), enolase, GAPDH, Mip, DnaJ | Deliberate safety exclusion; absolute one-way association (0/37 recovered antigens carry a human hit; p = 3.8 × 10⁻⁵) |
| Localisation list, Gram-positive `membrane` gap | 2 | CapB (*B. anthracis*), ActA (*L. monocytogenes*) | Configuration inconsistency, correctable |
| Upstream localisation misassignment (PSORTb) | 1 | LigA (*L. interrogans*) | Limit of the localisation predictor |
| Epitope prediction stage | 2 | ESAT-6 (*M. bovis*), L7/L12 (*Brucella*) | Prediction-level miss; toxicity hypothesis tested and rejected |

**VaxForge returns the same antigens on a virulent reference and three licensed vaccine strains.** Reverse vaccinology is ultimately judged on the genomes that deployed vaccines are built from, so we applied VaxForge at full proteome scale to three licensed live-attenuated brucellosis vaccine strains and set the result beside the virulent reference *B. suis* 1330: the cattle vaccines *B. abortus* S19 and RB51 and the small-ruminant vaccine *B. suis* S2, each on its own RefSeq proteome under default parameters (Table S2; Fig. S1). In every strain VaxForge recovered the outer-membrane and periplasmic protective antigens, Omp25, BP26/Omp28 and SodC throughout and Omp31 wherever the gene is present; *omp31* is deleted in *B. abortus*, so S19 and RB51 carry four ground-truth antigens rather than five (Cloeckaert et al., 2002). In every strain the single miss was the same protein, L7/L12, lost at the epitope stage exactly as in the virulent reference. The S2 vaccine, which shares the five-antigen set of *B. suis* 1330, reproduced its result closely (4 of 5; best rank 24 of 740; Fisher exact p = 0.015), and the two *B. abortus* vaccines reached 3 of 4 (best ranks 87 and 89; p = 0.055). The pipeline therefore returns the same answer on a virulent strain and its attenuated derivatives, and its one consistent omission is reproducible rather than stochastic.

**Host choice reshapes the predicted epitope repertoire in every virus.** VaxForge draws binding predictions from curated allele panels for human, mouse, bovine and pig, taken from IMGT/HLA and IPD-MHC (Table 5). For human, mouse and bovine, both classes run on real NetMHCpan and NetMHCIIpan models, including the BoLA-DRB3 class II panel for cattle; pig class I is real while its class II rests on a single validated allele. The panels are representative subsets rather than complete databases, and invalid allele names are dropped at runtime, so the reported epitopes are conservative.

**Table 5.** Host allele panels used for epitope prediction (representative subsets; full panels in `config/hosts.yaml`).

| Host | MHC-I panel (n) | MHC-II panel (n) | Predictor status |
|---|---|---|---|
| Human (HLA) | HLA-A and HLA-B (9) | HLA-DRB1 (6) | NetMHCpan and NetMHCIIpan, real |
| Mouse (H-2) | H-2 K/D (5) | H-2 IA/IE (5) | real |
| Bovine (BoLA) | BoLA class I (8) | BoLA-DRB3 (8) | real |
| Pig (SLA) | SLA class I (6) | SLA-DRB1 (1) | class I real; class II single allele |

The viral arm was built to exercise this layer across ten viruses that all four hosts naturally infect, spanning seven families, with the protein carrying the protective, patented or licensed target as ground truth (`PATENTED_VACCINES.md`). Viral envelope and capsid antigens are the intended vaccine targets yet are membrane-embedded, so a transmembrane-count filter tuned for bacteria discards them; we therefore made that filter soft for the viral profile and added a surface-exposure signal (Methods). With this in place VaxForge placed the protective antigen within the top quartile of candidates for eight of the ten viruses, first for Getah, Rift Valley fever and Japanese encephalitis and second for rabies, vesicular stomatitis and Venezuelan equine encephalitis (Table 6; Fig. 5). Against Vaxign-ML the ranking head-to-head was balanced, with four wins each and two ties, which is the honest outcome for two general rankers on viral proteomes. The one clear weakness is reovirus σ1 at rank 9 of 11: σ1 is a fibre protein with no membrane anchor and so receives no surface-exposure credit, whereas Vaxign-ML ranked it first, which is a real limitation of anchoring the viral signal on membrane topology.

The capability that separates VaxForge from both comparators is not the ranking but the host dimension, and it held across all ten viruses. For every virus the epitope set of the ground-truth antigen differed between host species, with pairwise host Jaccard overlap below one in every case and, for several viruses, host pairs sharing no epitopes at all (Table 6; Fig. 6; the rabies glycoprotein pairwise detail is in Fig. S2). NERVE 2.0 cannot process viral proteomes, and Vaxign-ML returns a single host-free protegenicity score, ranking rabies glycoprotein G fourth of five with no host resolution, so this per-host repertoire is unique to VaxForge among the three. Because each virus contributes a single protective antigen, the rank statistic is descriptive rather than inferential, and population coverage was reported as unavailable for non-human hosts rather than imputed.

**Table 6.** Host-specific epitope resolution and VaxForge–Vaxign-ML head-to-head across ten four-host viruses. Rank is the ground-truth antigen's position in each tool's candidate list; NERVE 2.0 is omitted because it cannot process viral proteomes. Jaccard is the mean (range) host-pair overlap of the ground-truth-antigen epitopes; below one means host choice changes the repertoire.

| Virus (species) | Ground-truth antigen | VaxForge rank | Vaxign-ML rank | Host-pair epitope Jaccard, mean (range) |
|---|---|---|---|---|
| Rabies virus (*Lyssavirus rabies*) | glycoprotein G | 2/5 | 4/5 | 0.62 (0.50–0.75) |
| Vesicular stomatitis Indiana virus (*Vesiculovirus indiana*) | glycoprotein G | 2/5 | 3/5 | 0.58 (0.33–0.78) |
| Encephalomyocarditis virus (*Cardiovirus A*) | capsid VP1 | 3/11 | 1/11 | 0.44 (0.00–0.78) |
| Rift Valley fever virus (*Phlebovirus riftense*) | glycoprotein Gn/Gc | 1/4 | 1/4 | 0.33 (0.00–1.00) |
| Venezuelan equine encephalitis virus (*Alphavirus VEE*) | E2 glycoprotein | 2/9 | 5/10 | 0.37 (0.00–1.00) |
| Getah virus (*Alphavirus getah*) | structural polyprotein (E2) | 1/3 | 1/3 | 0.50 (0.00–1.00) |
| Japanese encephalitis virus (*Orthoflavivirus japonicum*) | envelope protein E | 1/12 | 3/15 | 0.50 (0.00–1.00) |
| Mammalian orthoreovirus (*Orthoreovirus mammalis*) | σ1 attachment protein | 9/11 | 1/11 | 0.33 (0.00–1.00) |
| Vaccinia virus (*Orthopoxvirus vaccinia*) | L1R (IMV membrane) | 5/193 | 3/223 | 0.50 (0.00–1.00) |
| Cowpox virus (*Orthopoxvirus cowpox*) | L1R ortholog (CPXV099) | 8/207 | 3/233 | 0.50 (0.00–1.00) |

![Figure 5](figures/fig6_virus_antigen_rank.png)

**Figure 5.** Recovery of the protective or patented antigen across ten viruses, VaxForge versus Vaxign-ML. Each pair of bars is one virus; bar length is the ground-truth antigen's rank expressed as a percentile of that tool's candidate list (shorter is better), with the absolute rank labelled. The dotted line marks the top quartile. NERVE 2.0 is absent because it cannot process viral proteomes. Neither tool dominates: VaxForge leads on four viruses, Vaxign-ML on four, and the two tie on two.

![Figure 6](figures/fig7_virus_multihost.png)

**Figure 6.** Host selection reshapes the predicted epitope repertoire in every virus. Each row is one virus; the bar spans the minimum-to-maximum pairwise host Jaccard overlap of the ground-truth-antigen epitopes and the dot is the mean. A value of one would mean two hosts share an identical epitope set. Every virus falls below one, and several reach zero for at least one host pair, so the epitopes predicted for one host do not transfer to another. This host dimension is unavailable in both comparators.

---

**Discussion**

Benchmarks of reverse vaccinology tools have usually reported enrichment or classifier accuracy on curated positive and negative sets (Goodswen et al., 2023). We report instead how many independently validated protective antigens each tool returns from a whole proteome, and then, for the tool we built, where each missing antigen was lost. The second step is what turns a recall number into a usable statement. VaxForge recovered 37 of 48 antigens against 47 of 48 for Vaxign-ML, and on the headline number alone a reader would conclude that a 2020 descriptor-based classifier outperforms the pipeline. The attribution shows that six of the eleven misses were removed by a safety filter that Vaxign-ML does not apply at all, one by an upstream localisation error, two by a correctable configuration gap and two at the MHC scoring threshold. What remains after accounting for the deliberate exclusions is a far smaller difference in discrimination than the raw totals imply, and the direction of the host-homology association, no human homologue in any of the 37 recovered antigens, indicates a filter that is partitioning the ground-truth set on a biological axis rather than failing at random.

That filter is the central design choice of the pipeline and we defend it rather than apologise for it. All six antigens it removed are conserved housekeeping proteins, chaperones, glycolytic enzymes and an isomerase, whose protective activity in the literature is real but which necessarily carry human orthologues. Bacterial HSP60 cross-reacting with its human counterpart is an established route to autoimmune pathology (Matsuura et al., 2009), and shared peptides between pathogen and human heat shock proteins fall inside experimentally validated immunoreactive epitopes (Lucchese & Flöel, 2020), so excluding this class by default is a defensible prior for a pipeline intended to produce a candidate immunogen rather than a ranked list. Three honest qualifications follow. The 35% identity cut-off is a parameter choice, not a universal threshold. Sequence homology is not pathology, and we do not claim these six antigens would in fact cause autoimmunity, only that a cautious pipeline declines them. And the exclusion is a liability wherever a moonlighting antigen is the intended target, as in *S. suis*, where four of five ground-truth antigens are of this class and VaxForge consequently recovered one; a user pursuing such a target should relax the filter, and the pipeline exposes it as a single configuration value for exactly that reason.

Our enrichment finding converges with the comparator literature from the opposite direction. NERVE 2.0 reports fold enrichment near 16.6 for its best configuration and found Vaxign-ML to perform poorly on that metric, at 1.58 (Conte et al., 2024). We observed the mirror image on validated antigens: NERVE 2.0 reached fold enrichments of 25 to 274 while returning no protective antigen at all for four of nine species, and Vaxign-ML had the lowest enrichment and the highest recall. Both studies are measuring the same trade-off from different ends, which argues that neither metric should be reported alone. The point is sharpened by a design detail of the NERVE 2.0 benchmark: its ground-truth set deliberately excluded antigens annotated as cytosolic or of unknown localisation (Conte et al., 2024). Our set does not, so it contains precisely the moonlighting and conserved-housekeeping antigens that localisation-gated and homology-gated funnels lose, and our recall figures are correspondingly conservative relative to benchmarks built that way.

The comparison also needs placing against the current generation of antigen classifiers, which has moved from curated descriptors to learned representations. Vaxign-DL adds protein language model features to the Vaxign line (Chen et al., 2024), a descriptor-free ESM-2 pipeline outperformed descriptor-based classification on 9 of 10 bacterial species (Podda et al., 2025), and PLGDL combines sequence embeddings with geometric deep learning over predicted structures and confirmed a previously uncharacterised orthopoxvirus antigen in vivo (Zai et al., 2025). We did not benchmark against these, and the honest reason is that they answer a different question. Each returns a protegenicity score per protein; none returns host-specific epitopes, a construct or an immunisation readout, and PLGDL additionally requires predicted structures and substantial compute. VaxForge is not competing on classification accuracy, where a modern language model will plausibly beat its funnel, and a productive combination is visible: a learned classifier is the natural replacement for the scoring stage of the funnel, feeding the host-resolved immunological layer that the classifiers themselves do not provide. Treating the two as alternatives rather than as stages is, in our view, what has kept the field's attention on the ranking step.

The host dimension is the contribution we consider least substitutable. NetMHCpan and NetMHCIIpan already support many non-human alleles (Reynisson et al., 2020), yet reverse vaccinology pipelines apply them against a single human reference, and NERVE 2.0's epitope module is explicitly built on human HLA supertypes (Conte et al., 2024). Because all ten viruses infect the four hosts examined (Cox et al., 1977; Fooks et al., 2017), we could compute a per-host repertoire for each, and the epitope sets of the ground-truth antigens differed between species in every case, with some host pairs sharing nothing. For the veterinary and zoonotic targets that carry much of the practical vaccine burden (Fooks et al., 2017; Hampson et al., 2015), a human-only epitope layer does not merely lose precision, it reports a repertoire for the wrong species.

Several limitations bound these claims, beyond the three qualifications already given for the homology filter. The Gram-positive localisation gap is a defect, not a trade-off: the permitted-localisation list encodes Gram-negative architecture and cost two legitimate surface antigens, and correcting it requires Gram-aware lists and a full re-run of all ten species, so the recall reported here should be read as a floor. One miss follows an upstream PSORTb misassignment that we can report but not repair. Two misses occur at MHC scoring thresholds that we did not tune. The default toxicity filter can discard peptides from toxin antigens such as listeriolysin O, a conservative choice left untuned for the benchmark; most B-cell epitopes are discontinuous and a linear-epitope pipeline does not model them (Michalik et al., 2016). Area-under-curve is reported only as a secondary measure with an explicit small-sample caveat, because four to five positives per pathogen cannot resolve AUC differences. Whole-proteome VaxForge runs remain impractical on CPU-only hardware, so the meningococcal comparison used matched subsamples, fair across tools, while full proteomes were run for the faster comparators. Finally, and most importantly, nothing here is experimental validation: every result is a prediction about predictions, and the ground-truth sets themselves reflect protection observed in particular animal models and immunisation protocols.

Three directions follow directly. Correcting the Gram-aware localisation lists and re-running the panel would separate the correctable component of the recall gap from the deliberate one. Substituting a protein language model classifier for the scoring stage of the funnel would test whether the end-to-end architecture and modern discrimination are complementary, as we expect. And the cross-species divergence of the viral glycoprotein epitopes is a concrete, testable hypothesis: the predicted repertoires for cattle and human rabies glycoprotein G share two thirds of their epitopes at most, which binding and immunisation experiments could confirm or refute.

---

**Methods**

**Candidate discovery and antigen funnel.** VaxForge is deterministic and reads every cutoff from a single configuration file, so that no threshold is hard-coded (`config/thresholds.yaml`); input type, whether a proteome, coding sequences or a genome, was detected automatically. Proteins were aligned to the curated virulence factor database VFDB (Liu et al., 2022) with DIAMOND blastp (Buchfink et al., 2021) at an e-value of at most 1 × 10⁻⁵, sequence identity of at least 30% and coverage of at least 50%; this evidence contributed to the candidacy score and was not used as a gate. Subcellular localisation was predicted with PSORTb 3.0 (Yu et al., 2010), with DeepLoc 2.0 (Thumuluri et al., 2022) as a fallback, retaining extracellular, outer-membrane, cell-wall and periplasmic proteins and removing cytoplasmic ones (hard filter). Transmembrane topology was predicted with TMHMM-2.0 (Krogh et al., 2001) and proteins with more than one helix were removed (hard filter). Human homology was assessed with DIAMOND blastp against the human proteome and proteins above 35% identity were excluded for autoimmune safety (hard filter). Signal peptides (SignalP-5.0; Almagro Armenteros et al., 2019, probability ≥ 0.5), antigenicity (IApred; Miles et al., 2025), cross-strain conservation (≥ 80% when multiple strains were supplied) and adhesin probability (≥ 0.51) entered as soft, score-weighted terms. Surviving proteins were ranked by a weighted candidacy score.

**Host-specific epitope prediction.** For the selected host, class I binding was predicted with NetMHCpan-4.1 (Reynisson et al., 2020) for peptides of 8 to 11 residues, with strong and weak binders defined at percentile ranks of ≤ 0.5% and ≤ 2%. Antigen processing was added from NetCTL-1.2 as allele-independent proteasomal cleavage and TAP transport (Larsen et al., 2005), and class I immunogenicity from the IEDB model (Calis et al., 2013). Class II binding was predicted with NetMHCIIpan (Reynisson et al., 2020) for 15-mers with a 9-residue core, at percentile ranks of ≤ 2% and ≤ 10%, and helper epitopes were further screened for interferon-γ induction with IFNepitope2 (score ≥ 0.49; Dhall et al., 2024) in human and mouse. Linear B-cell epitopes were predicted with BepiPred-3.0 (Clifford et al., 2022).

**Epitope safety and population coverage.** Retained epitopes were required to be non-allergenic, judged by a six-residue identity match to known allergens (FAO/WHO, 2001), and non-toxic, judged by ToxinPred2 with peptides above a machine-learning score of 0.6 removed (Sharma et al., 2022), both as hard filters; epitope self-mimicry to the host proteome, defined as an exact 9-residue match, was flagged softly. Population coverage of the selected epitope set across human HLA alleles was computed with the IEDB population-coverage tool (Bui et al., 2006), with a target of at least 90%.

**Multi-epitope construct and in silico immunisation.** The retained epitopes were joined into a single construct with AAY linkers between cytotoxic epitopes and KK linkers around B-cell epitopes (Schubert & Kohlbacher, 2016), GPGPG linkers between helper epitopes (Livingston et al., 2002) and an N-terminal β-defensin adjuvant (Biragyn et al., 2002). Physicochemical properties of the construct were computed with ProtParam (Gasteiger et al., 2005), and the construct was passed to an in silico immunisation module that simulates the antibody and isotype response over time.

**Benchmark pathogens and ground-truth antigens.** The bacterial arm used serogroup B *N. meningitidis* MC58 (UniProt proteome UP000000425, 2,003 proteins; RefSeq NC_003112.2) with the four 4CMenB antigens fHbp (Q9JXV4), NadA (Q9JXK7), NHBA (Q7DD37) and PorA (P0DH58) as ground truth, together with nine further species at full proteome scale. The first three were *S. aureus* NCTC 8325 (UP000008816, 2,889 proteins), *L. monocytogenes* EGD-e (UP000000817, 2,844) and *S.* Typhimurium LT2 (UP000001014, 4,533), each with five protective or candidate antigens. The remaining six were zoonotic bacteria that all four mammalian hosts can infect, selected to span genome size and cell biology: *Coxiella burnetii* RSA 493 (UP000002671, 1,812 proteins; obligate intracellular), *Streptococcus suis* serotype 2 BM407 (UP000009077, 1,922), *Brucella suis* 1330 (UP000007104, 3,249), *Mycobacterium bovis* AF2122/97 (UP000001419, 4,088), *Leptospira interrogans* serovar Lai 56601 (UP000001408, 3,676; a spirochaete) and *Bacillus anthracis* Ames Ancestor (UP000000594, 5,493 proteins, chromosome plus the pXO1 and pXO2 plasmids that carry the toxin and capsule antigens), with four or five protective antigens each. The protective antigens of all nine species, with their evidence and DOIs, are listed in `PROTECTIVE_ANTIGENS.md` and their accessions in `PATHOGEN_ACCESSIONS.md`. The viral arm used ten viruses that all four mammalian hosts naturally infect, across seven families: rabies virus (*Lyssavirus rabies*; NC_001542, glycoprotein G at NP_056796) and vesicular stomatitis Indiana virus (*Vesiculovirus indiana*; NC_001560, glycoprotein G) of the Rhabdoviridae; encephalomyocarditis virus (*Cardiovirus A*; NC_001479, capsid VP1) of the Picornaviridae; Venezuelan equine encephalitis virus (NC_001449, E2) and Getah virus (NC_006558, structural polyprotein carrying E2) of the Togaviridae; Japanese encephalitis virus (NC_001437, envelope E) of the Flaviviridae; Rift Valley fever virus (NC_014397, glycoprotein Gn/Gc) of the Phenuiviridae; mammalian orthoreovirus (σ1 attachment protein) of the Reoviridae; and vaccinia (NC_006998, L1R) and cowpox (NC_003663, L1R ortholog) viruses of the Poxviridae. Each ground-truth antigen carries a protective, licensed or patented vaccine target, with patents and evidence listed in `PATENTED_VACCINES.md`. For the polyprotein and segmented viruses, mature peptides were extracted where annotated, and for the two orthopoxviruses a single RefSeq reference was used with the cowpox L1R identified by orthology to the vaccinia L1R. All accessions were verified against NCBI and UniProt (`PATHOGEN_ACCESSIONS.md`). For meningococcus a fixed-seed subsample of the proteome (64 and 204 proteins, four protective antigens plus background) was used, because whole-proteome VaxForge runs are impractical on CPU-only hardware, and the same subsample was given to every tool.

**Comparator tools and their roles.** VaxForge was compared with the two reverse vaccinology tools that are current, locally installable and most widely used for bacterial antigen prediction, so that the comparison reflects what a user would actually run (Table 1). Vaxign-ML (Ong et al., 2020) is a supervised machine-learning classifier that takes a protein proteome and returns a protegenicity probability for each protein, ranking candidates but producing no epitopes, host information or construct. NERVE 2.0 (Conte et al., 2024) is an end-to-end bacterial reverse vaccinology environment that takes a bacterial proteome and a Gram label and returns a shortlist of candidate proteins after a cascade of hard selection filters, followed by an optional epitope-prediction module that reports linear MHC-I/II binders and population coverage for human HLA supertypes. VaxForge takes a proteome, coding sequences or a genome, for bacteria or viruses, and returns ranked antigens together with host-specific epitopes, a multi-epitope construct and an in silico immunisation readout. Both comparators were installed locally in containers and run with default parameters on identical input, with no threshold tuning; recall over ground-truth antigens and the antigen ranks were read from each tool's native output. Protein language model classifiers (Chen et al., 2024; Podda et al., 2025; Zai et al., 2025) were not included because they return a protegenicity score only and, for the structure-based framework, require predicted structures and GPU compute; the rationale is set out in the Discussion.

**Attribution of missed antigens.** For each ground-truth antigen that VaxForge did not return, the pipeline's intermediate outputs were inspected to identify the first step at which it was removed, using the recorded PSORTb localisation call and score, the TMHMM helix count, and DIAMOND blastp identity against the human Swiss-Prot proteome (`tools/db/human.dmnd`, e-value 1 × 10⁻⁴). Antigens that passed all three hard filters but produced no candidate peptide were assigned to the epitope stage. Association between host-homology status and recovery was tested with the Fisher exact test on the two-by-two table of recovered against human-homologous antigens across all 48 ground-truth antigens. For the two epitope-stage misses the toxicity hypothesis was tested by comparing the fraction of flagged 9-mers in the missed antigens against a recovered control of comparable length from the same species. Per-antigen values are in `results/miss_attribution.json` and the full reasoning in `MISS_ATTRIBUTION.md`.

**Statistical analysis.** Enrichment of ground-truth antigens among candidates was tested with the Fisher exact test, and antigen ranks with a distribution-free permutation test. Discrimination was summarised with ROC-AUC and bootstrap confidence intervals, and AUC differences with the DeLong test; because each pathogen carries only four to five positives, AUC is reported as secondary with a small-sample caveat. Overlap of host epitope sets was quantified with the Jaccard index. All tests are implemented in `benchmark/stats.py`.

**Pipeline scale-safety validation.** During the full-proteome bacterial runs, three scale-dependent defects that had artefactually lowered VaxForge recall were identified and corrected. The localisation step originally submitted every protein in a single container call under a fixed timeout, which for the largest proteome returned empty and forced a heuristic fallback that misassigned surface antigens; it now runs in sub-batches of at most 800 proteins. The class I step originally called NetMHCpan once with the concatenated peptides of all proteins (about 1.37 million unique 9-mers at whole-proteome scale), which collapsed to near-zero coverage; it now processes peptides in chunks of 20,000 in parallel across cores. A third defect appeared only at the largest proteomes: the parallel class I workers shared one scratch directory, so concurrent NetMHCpan instances could deadlock on each other and stall a run indefinitely, which was resolved by giving each worker an isolated temporary directory and a bounded per-chunk timeout. Control experiments confirmed that the toxicity and allergenicity filters were not responsible, and the affected species were re-run after the fixes. The miss attribution reported here was performed after these corrections and found no residual scale-dependent cause. Full detail is in `MANUSCRIPT_METHODS.md`, section B.

**Data and code availability.** Benchmark scripts, results and figures are under `benchmark/`; pipeline method citations are in `vaxforge/citations.py`; all thresholds are in `config/thresholds.yaml`.

---

**Supplementary material**

**Table S1.** The ground-truth protective antigens, by species, with accession, licensure status and the verified evidence for protective or candidate status (full evidence and DOIs in `PROTECTIVE_ANTIGENS.md`). Status: *licensed* = component of an approved human vaccine; *vet-licensed strain* = genus has approved live/whole-cell veterinary vaccines while the subunit here is a candidate; *clinical/preclinical candidate* = experimental protective or immunodominant antigen.

| Species | Antigen (gene) | Accession | Status | VaxForge | Vaxign-ML | NERVE 2.0 | Evidence (DOI) |
|---|---|---|---|---|---|---|---|
| *N. meningitidis* B | fHbp | Q9JXV4 | Licensed, 4CMenB (Bexsero) | #2 | #7 | filt. | Giuliani 2006, 10.1073/pnas.0603940103 |
|  | NadA | Q9JXK7 | Licensed, 4CMenB | #4 | #1 | filt. | Giuliani 2006, 10.1073/pnas.0603940103 |
|  | NHBA | Q7DD37 | Licensed, 4CMenB | #3 | #2 | filt. | Giuliani 2006, 10.1073/pnas.0603940103 |
|  | PorA | P0DH58 | Licensed, 4CMenB (OMV) | #1 | #4 | #2 | Pizza 2000, 10.1126/science.287.5459.1816 |
| *S. aureus* | IsdB | Q2FZF0 | Clinical candidate (V710) | #4 | #46 | filt. | Kuklin 2006, 10.1128/iai.74.4.2215-2223.2006 |
|  | IsdA | Q2FZE9 | Preclinical candidate | #10 | #387 | filt. | Kim 2010, 10.1016/j.vaccine.2010.02.097 |
|  | ClfA | Q2G015 | Preclinical candidate | #23 | #248 | #3 | Josefsson 2001, 10.1086/324430 |
|  | Hla | Q2G1X0 | Preclinical candidate | #9 | #29 | #24 | Adhikari 2012, 10.1371/journal.pone.0038567 |
|  | SpA | P02976 | Preclinical candidate | #8 | #178 | #28 | Kim 2010, 10.1084/jem.20092514 |
| *L. monocytogenes* | LLO | P13128 | Preclinical candidate | #7 | #19 | filt. | Phelps 2020, 10.1016/j.vaccine.2020.06.049 |
|  | p60 | P21171 | Preclinical (CD8) | #1 | #70 | #30 | Harty 1995, 10.4049/jimmunol.154.9.4642 |
|  | InlA | P0DJM0 | Virulence factor / candidate | #2 | #214 | filt. | Gaillard 1991, 10.1016/0092-8674(91)90009-n |
|  | InlB | P0DQD2 | Virulence factor / candidate | #15 | #143 | #72 | Dramsi 1995, 10.1111/j.1365-2958.1995.tb02297.x |
|  | ActA | P33379 | Passive-protection candidate | ✗ | #210 | #38 | Asano 2016, 10.1038/srep39628 |
| *S.* Typhimurium | OmpD | P37592 | Preclinical candidate | #63 | #28 | filt. | Gil-Cruz 2009, 10.1073/pnas.0812431106 |
|  | OmpC | P0A263 | Preclinical candidate | #51 | #30 | #70 | Secundino 2006, 10.1111/j.1365-2567.2005.02263.x |
|  | OmpA | P02936 | Preclinical candidate | #20 | #185 | #51 | Domínguez-Medina 2020, 10.1038/s41467-020-14655-9 |
|  | FliC | P06179 | Preclinical candidate | #8 | #127 | #6 | Bumann 2014, 10.3389/fimmu.2014.00381 |
|  | SseB | Q7BVH7 | Preclinical candidate | #12 | #897 | filt. | Lee 2017, 10.4049/jimmunol.1601357 |
| *C. burnetii* | Com1 | H7C7D7 | Preclinical (adoptive transfer) | #42 | #64 | filt. | Wang 2011, 10.1186/1471-2172-12-52 |
|  | Mip | P51752 | Immunodominant candidate | ✗ | #28 | filt. | Wang 2011, 10.1186/1471-2172-12-52 |
|  | OmpH | Q83DT1 | Immunodominant candidate | #177 | #713 | filt. | Wang 2011, 10.1186/1471-2172-12-52 |
|  | GroEL | P19421 | Immunodominant (weak) | ✗ | #8 | filt. | Wang 2011, 10.1186/1471-2172-12-52 |
| *S. suis* | Suilysin | A0A0H3MU02 | Preclinical candidate | #6 | #52 | filt. | Du 2013, 10.1016/j.vaccine.2013.07.004 |
|  | Enolase | A0A0H3N5B2 | Preclinical (full protection) | ✗ | #150 | filt. | Feng 2009, 10.1086/644602 |
|  | 6PGD | A0A0H3N193 | Preclinical candidate | ✗ | #68 | filt. | Zhang 2023, 10.3390/vetsci10060383 |
|  | GAPDH | A0A0H3MSV3 | Moonlighting candidate | ✗ | #220 | filt. | Feng 2009, 10.1086/644602 |
|  | DnaJ | A0A0H3MXM0 | Candidate | ✗ | #482 | filt. | Zhang 2023, 10.3390/vetsci10060383 |
| *B. suis* | Omp25 | Q45689 | Vet-licensed strain / candidate | #102 | #101 | #3 | Cloeckaert 2002, 10.1016/s0378-1135(02)00211-0 |
|  | Omp31 | P0A3U5 | Preclinical (protective) | #36 | #218 | #4 | Cloeckaert 2002, 10.1016/s0378-1135(02)00211-0 |
|  | BP26 | P0A3U9 | Immunodominant / diagnostic | #181 | #255 | filt. | Kim 2013, 10.1016/j.jmb.2013.01.015 |
|  | SodC | P66827 | Preclinical (protective) | #26 | #589 | filt. | Oñate 2003, 10.1128/IAI.71.9.4857-4861.2003 |
|  | L7/L12 | P0A469 | Preclinical (protective) | ✗ | #430 | filt. | Luo 2006, 10.1128/IAI.74.5.2734-2741.2006 |
| *M. bovis* | ESAT-6 | P0A565 | Clinical candidate (BCG-related) | ✗ | #53 | filt. | Brandt 2000, 10.1128/IAI.68.2.791-795.2000 |
|  | Ag85B | P0C2T2 | Clinical candidate (H1/H56) | #1 | #38 | filt. | Davila 2010, 10.1155/2010/208371 |
|  | CFP-10 | P0A567 | Candidate / diagnostic | #27 | #70 | filt. | Davila 2010, 10.1155/2010/208371 |
|  | MPB64 | P0A5Q5 | Immunodominant / diagnostic | #202 | #13 | filt. | Brandt 2000, 10.1128/IAI.68.2.791-795.2000 |
|  | HspX | P0A5B8 | Immunodominant candidate | #114 | #342 | filt. | Davila 2010, 10.1155/2010/208371 |
| *L. interrogans* | LigA | Q8EYU4 | Vet-licensed strain / candidate | ✗ | #227 | filt. | Faisal 2007, 10.1016/j.vaccine.2007.10.029 |
|  | LigB | Q8EZS3 | Preclinical candidate | #171 | #18 | filt. | Faisal 2007, 10.1016/j.vaccine.2007.10.029 |
|  | LipL32 | O34094 | Preclinical (protective) | #145 | #788 | filt. | Feng 2009, 10.1590/s0100-879x2009005000013 |
|  | LipL41 | Q8F8E1 | Preclinical (protective) | #107 | #83 | filt. | Feng 2009, 10.1590/s0100-879x2009005000013 |
|  | OmpL1 | G1UB30 | Preclinical (protective) | #150 | #588 | filt. | Feng 2009, 10.1590/s0100-879x2009005000013 |
| *B. anthracis* | PA (pagA) | P13423 | Licensed (BioThrax/AVA, AVP) | #2 | #13 | #3 | Rhie 2005, 10.1016/j.femsim.2005.05.009 |
|  | LF (lef) | P15917 | Preclinical (protective) | #9 | #108 | filt. | Baillie 2010, 10.1016/j.vaccine.2010.07.075 |
|  | EF (cya) | P40136 | Toxin component / candidate | #6 | #119 | filt. | Baillie 2010, 10.1016/j.vaccine.2010.07.075 |
|  | CapB (γ-PGA) | P19580 | Preclinical (protective adjunct) | ✗ | #1196 | filt. | Schneerson 2003, 10.1073/pnas.1633512100 |
|  | EA1 (S-layer) | P94217 | Immunodominant candidate | #65 | #2 | #2 | Zhang 2007, 10.1016/j.vaccine.2007.11.072 |

**Table S2.** VaxForge applied to the virulent *Brucella* reference and three licensed live-attenuated vaccine strains at full proteome scale (default parameters, identical pipeline). Each cell is the antigen's rank in the VaxForge candidate list; "–" marks an antigen absent from that genome (*omp31* is deleted in *B. abortus*) and ✗ an antigen not recovered. Recall counts recovered ground-truth antigens; the Fisher exact p tests their enrichment among candidates.

| Strain (status) | Species | Proteins | Omp25 | Omp31 | BP26 | SodC | L7/L12 | Recall | Fisher p |
|---|---|---|---|---|---|---|---|---|---|
| 1330 (virulent reference) | *B. suis* | 3,249 | #102 | #36 | #181 | #26 | ✗ | 4/5 | 0.026 |
| RB51 (licensed, cattle) | *B. abortus* | 2,971 | #156 | – | #87 | #524 | ✗ | 3/4 | 0.055 |
| S19 (licensed, cattle) | *B. abortus* | 2,964 | #159 | – | #89 | #526 | ✗ | 3/4 | 0.055 |
| S2 (licensed, small ruminant) | *B. suis* | 2,979 | #54 | #24 | #94 | #103 | ✗ | 4/5 | 0.015 |

![Figure S1](figures/figS1_brucella_vaccine.png)

**Figure S1.** VaxForge gives a consistent answer on the virulent *Brucella* reference and three licensed live-attenuated vaccine strains. Each row is one protective antigen and each point is one strain's rank for that antigen in the VaxForge candidate list, so a point further left is a better rank; the candidate-list size n for each strain is in the legend. A cross in the "not recovered" column marks an antigen the strain did not return. The outer-membrane and periplasmic antigens, Omp25, Omp31, BP26/Omp28 and SodC, are recovered in every strain that carries them; Omp31 is absent from the two *B. abortus* vaccines because the gene is deleted in that species (Cloeckaert et al., 2002); and the ribosomal antigen L7/L12 is lost at the epitope-prediction stage in all four strains, the same reproducible call seen in the virulent reference.

![Figure S2](figures/fig5_multihost_jaccard.png)

**Figure S2.** Host selection reshapes the rabies glycoprotein G epitope repertoire, shown pairwise. Each bar is one pair of host species (for example human–pig), and its value is the Jaccard overlap between the glycoprotein G epitopes predicted for those two hosts: 1 would mean the two hosts share an identical epitope set and 0 means they share none. Every pair falls below one, so the epitopes predicted for one host do not fully transfer to another and the choice of host changes the predicted repertoire.

---

**References**

Almagro Armenteros, J. J., Tsirigos, K. D., Sønderby, C. K., Petersen, T. N., Winther, O., Brunak, S., von Heijne, G., & Nielsen, H. (2019). SignalP 5.0 improves signal peptide predictions using deep neural networks. *Nature Biotechnology, 37*, 420–423. https://doi.org/10.1038/s41587-019-0036-z

Biragyn, A., Ruffini, P. A., Leifer, C. A., Klyushnenkova, E., Shakhov, A., Chertov, O., Coscia, M., Boyd, L. F., Xu, Y., & Kwak, L. W. (2002). Toll-like receptor 4-dependent activation of dendritic cells by β-defensin 2. *Science, 298*(5595), 1025–1029. https://doi.org/10.1126/science.1075565

Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods, 18*, 366–368. https://doi.org/10.1038/s41592-021-01101-x

Bui, H. H., Sidney, J., Dinh, K., Southwood, S., Newman, M. J., & Sette, A. (2006). Predicting population coverage of T-cell epitope-based diagnostics and vaccines. *BMC Bioinformatics, 7*, 153. https://doi.org/10.1186/1471-2105-7-153

Calis, J. J. A., Maybeno, M., Greenbaum, J. A., Weiskopf, D., De Silva, A. D., Sette, A., Keşmir, C., & Peters, B. (2013). Properties of MHC class I presented peptides that enhance immunogenicity. *PLoS Computational Biology, 9*(10), e1003266. https://doi.org/10.1371/journal.pcbi.1003266

Chen, Y., Zhang, Y., & He, Y. (2024). Enhancing Vaxign-DL for vaccine candidate prediction with added ESM-generated features. *bioRxiv* [Preprint]. https://doi.org/10.1101/2024.09.04.611295

Clifford, J. N., Høie, M. H., Deleuran, S., Peters, B., Nielsen, M., & Marcatili, P. (2022). BepiPred-3.0: Improved B-cell epitope prediction using protein language models. *Protein Science, 31*(12), e4497. https://doi.org/10.1002/pro.4497

Cloeckaert, A., Vizcaíno, N., Paquet, J. Y., Bowden, R. A., & Elzer, P. H. (2002). Major outer membrane proteins of *Brucella* spp.: past, present and future. *Veterinary Microbiology, 90*(1–4), 229–247. https://doi.org/10.1016/s0378-1135(02)00211-0

Conte, A., Gulmini, N., Costa, F., Cartura, M., Bröhl, F., Patanè, F., & Filippini, F. (2024). NERVE 2.0: boosting the new enhanced reverse vaccinology environment via artificial intelligence and a user-friendly web interface. *BMC Bioinformatics, 25*, 378. https://doi.org/10.1186/s12859-024-06004-0

Cox, J. H., Dietzschold, B., & Schneider, L. G. (1977). Rabies virus glycoprotein. II. Biological and serological characterization. *Infection and Immunity, 16*(3), 754–759. https://doi.org/10.1128/iai.16.3.754-759.1977

Dhall, A., Patyal, S., & Raghava, G. P. S. (2024). A hybrid method for discovering interferon-gamma inducing peptides in human and mouse. *Scientific Reports, 14*(1), 26859. https://doi.org/10.1038/s41598-024-77957-8

FAO/WHO. (2001). *Evaluation of allergenicity of genetically modified foods. Report of a Joint FAO/WHO Expert Consultation.* https://www.fao.org/3/y0820e/y0820e00.htm

Fooks, A. R., Cliquet, F., Finke, S., Freuling, C., Hemachudha, T., Mani, R. S., Müller, T., Nadin-Davis, S., Picard-Meyer, E., Wilde, H., & Banyard, A. C. (2017). Rabies. *Nature Reviews Disease Primers, 3*, 17091. https://doi.org/10.1038/nrdp.2017.91

Gasteiger, E., Hoogland, C., Gattiker, A., Duvaud, S., Wilkins, M. R., Appel, R. D., & Bairoch, A. (2005). Protein identification and analysis tools on the ExPASy server. In J. M. Walker (Ed.), *The Proteomics Protocols Handbook* (pp. 571–607). Humana Press. https://doi.org/10.1385/1-59259-890-0:571

Giuliani, M. M., Adu-Bobie, J., Comanducci, M., Aricò, B., Savino, S., Santini, L., Brunelli, B., Bambini, S., Biolchi, A., Capecchi, B., Cartocci, E., Ciucchi, L., Di Marcello, F., Ferlicca, F., Galli, B., Luzzi, E., Masignani, V., Serruto, D., Veggi, D., … Pizza, M. (2006). A universal vaccine for serogroup B meningococcus. *Proceedings of the National Academy of Sciences, 103*(29), 10834–10839. https://doi.org/10.1073/pnas.0603940103

Goodswen, S. J., Kennedy, P. J., & Ellis, J. T. (2023). A guide to current methodology and usage of reverse vaccinology towards in silico vaccine discovery. *FEMS Microbiology Reviews, 47*(2), fuad004. https://doi.org/10.1093/femsre/fuad004

Hampson, K., Coudeville, L., Lembo, T., Sambo, M., Kieffer, A., Attlan, M., Barrat, J., Blanton, J. D., Briggs, D. J., Cleaveland, S., Costa, P., Freuling, C. M., Hiby, E., Knopf, L., & Leanes, F. (2015). Estimating the global burden of endemic canine rabies. *PLoS Neglected Tropical Diseases, 9*(4), e0003709. https://doi.org/10.1371/journal.pntd.0003709

Krogh, A., Larsson, B., von Heijne, G., & Sonnhammer, E. L. L. (2001). Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *Journal of Molecular Biology, 305*(3), 567–580. https://doi.org/10.1006/jmbi.2000.4315

Larsen, M. V., Lundegaard, C., Lund, O., & Nielsen, M. (2005). An integrative approach to CTL epitope prediction: a combined algorithm integrating MHC class I binding, TAP transport efficiency, and proteasomal cleavage predictions. *European Journal of Immunology, 35*(8), 2295–2303. https://doi.org/10.1002/eji.200425811

Liu, B., Zheng, D., Zhou, S., Chen, L., & Yang, J. (2022). VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Research, 50*(D1), D912–D917. https://doi.org/10.1093/nar/gkab1107

Livingston, B., Crimi, C., Newman, M., Higashimoto, Y., Appella, E., Sidney, J., & Sette, A. (2002). A rational strategy to design multiepitope immunogens based on multiple Th lymphocyte epitopes. *Journal of Immunology, 168*(11), 5499–5506. https://doi.org/10.4049/jimmunol.168.11.5499

Lucchese, G., & Flöel, A. (2020). SARS-CoV-2 and Guillain-Barré syndrome: molecular mimicry with human heat shock proteins as potential pathogenic mechanism. *Cell Stress and Chaperones, 25*(5), 731–735. https://doi.org/10.1007/s12192-020-01145-6

Matsuura, E., Kobayashi, K., Matsunami, Y., Shen, L., Quan, N., Makarova, M., Suchkov, S. V., Ayada, K., Oguma, K., & Lopez, L. R. (2009). Autoimmunity, infectious immunity, and atherosclerosis. *Journal of Clinical Immunology, 29*(6), 714–721. https://doi.org/10.1007/s10875-009-9333-5

Michalik, M., Djahanshiri, B., Leo, J. C., & Linke, D. (2016). Reverse vaccinology: the pathway from genomes and epitope predictions to tailored recombinant vaccines. In *Vaccine Design* (Methods in Molecular Biology, Vol. 1403, pp. 87–106). Springer. https://doi.org/10.1007/978-1-4939-3387-7_4

Miles, S., Menafra, G., Iriarte, A., & Chabalgoity, J. A. (2025). IApred: A versatile open-source tool for predicting protein antigenicity across diverse pathogens. *ImmunoInformatics*. https://doi.org/10.1016/j.immuno.2025.100014

Ong, E., Wang, H., Wong, M. U., Seetharaman, M., Valdez, N., & He, Y. (2020). Vaxign-ML: supervised machine learning reverse vaccinology model for improved prediction of bacterial protective antigens. *Bioinformatics, 36*(10), 3185–3191. https://doi.org/10.1093/bioinformatics/btaa119

Pizza, M., Scarlato, V., Masignani, V., Giuliani, M. M., Aricò, B., Comanducci, M., Jennings, G. T., Baldi, L., Bartolini, E., Capecchi, B., Galeotti, C. L., Luzzi, E., Manetti, R., Marchetti, E., Mora, M., Nuti, S., Ratti, G., Santini, L., Savino, S., … Rappuoli, R. (2000). Identification of vaccine candidates against serogroup B meningococcus by whole-genome sequencing. *Science, 287*(5459), 1816–1820. https://doi.org/10.1126/science.287.5459.1816

Podda, M., Savojardo, C., Martelli, P. L., Casadio, R., Sîrbu, A., Priami, C., & Brozzi, A. (2025). A descriptor-free machine learning framework to improve antigen discovery for bacterial pathogens. *PLOS ONE, 20*(6), e0323895. https://doi.org/10.1371/journal.pone.0323895

Reynisson, B., Alvarez, B., Paul, S., Peters, B., & Nielsen, M. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research, 48*(W1), W449–W454. https://doi.org/10.1093/nar/gkaa379

Schubert, B., & Kohlbacher, O. (2016). Designing string-of-beads vaccines with optimal spacers. *Genome Medicine, 8*, 9. https://doi.org/10.1186/s13073-016-0263-6

Sharma, N., Naorem, L. D., Jain, S., & Raghava, G. P. S. (2022). ToxinPred2: an improved method for predicting toxicity of proteins. *Briefings in Bioinformatics, 23*(5), bbac174. https://doi.org/10.1093/bib/bbac174

Thumuluri, V., Almagro Armenteros, J. J., Johansen, A. R., Nielsen, H., & Winther, O. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. *Nucleic Acids Research, 50*(W1), W228–W234. https://doi.org/10.1093/nar/gkac278

Yu, N. Y., Wagner, J. R., Laird, M. R., Melli, G., Rey, S., Lo, R., Dao, P., Sahinalp, S. C., Ester, M., Foster, L. J., & Brinkman, F. S. L. (2010). PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. *Bioinformatics, 26*(13), 1608–1615. https://doi.org/10.1093/bioinformatics/btq249

Zai, X., Zhao, Y., Wang, X., Leng, M., Lu, M., Yang, Y., Zhao, X., Li, R., Li, Y., Zhang, Y., Zhang, J., Li, D., Ren, H., Xu, J., & Chen, W. (2025). Integrating protein language and geometric deep learning models for enhanced vaccine antigen prediction. *Nature Communications, 17*(1), 1033. https://doi.org/10.1038/s41467-025-67778-2
