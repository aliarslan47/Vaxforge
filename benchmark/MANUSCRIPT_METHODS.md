# VaxForge — makale-hazır Methods bölümleri

Bu dosya, benchmark "yazım aşaması" maddelerini **doğrudan makaleye taşınabilir** düzyazı olarak
toplar. İki bölüm: **(A) Rationale for benchmark selection** (item 3, kaynak
`SELECTION_RATIONALE.md` + `PROTECTIVE_ANTIGENS.md`) ve **(B) Pipeline scale-safety validation**
(item 2, kaynak `BENCHMARK_SUMMARY.md` §metodoloji). Metin İngilizce (yayın dili); tüm atıflar
DOI-doğrulanmış (`PROTECTIVE_ANTIGENS.md`, 2026-09-20). İlke: uydurma yok ([[vaxforge-literature-rule]]).

---

## A. Rationale for benchmark selection

### A.1 Pathogen selection

We evaluated VaxForge against two established reverse-vaccinology tools, NERVE 2.0 (Conte et al.,
2024) and Vaxign-ML (Ong et al., 2020), across a benchmark deliberately spanning both a bacterial
and a viral arm. Serogroup B *Neisseria meningitidis* (MenB) anchors the bacterial arm because it
is the pathogen on which reverse vaccinology was founded: Pizza et al. (2000) derived vaccine
candidates directly from the whole-genome sequence, defining the paradigm. MenB therefore provides
the field-standard test bed, and its licensed four-component vaccine (4CMenB) supplies an
experimentally validated, multi-antigen ground truth — recombinant fHbp, NadA and NHBA together
with a PorA-containing outer-membrane vesicle (Giuliani et al., 2006) — that is ideal for a
**recall@antigen** metric (how many known protective antigens a tool recovers).

To move beyond a single exemplar, the bacterial arm was extended to three additional zoonotic
species chosen along three orthogonal axes: increasing genome size (*Staphylococcus aureus*
NCTC 8325, 2.8 Mb → *Listeria monocytogenes* EGD-e, 2.9 Mb → *Salmonella enterica* serovar
Typhimurium LT2, 4.9 Mb), Gram balance (two Gram-positive, one Gram-negative, exercising the
Gram-specific localization path of PSORTb), and the availability of literature-established
protective/candidate antigens for a measurable recall metric.

The viral arm is built on rabies virus (RABV), whose single surface glycoprotein G is the only
structural protein that induces virus-neutralizing antibodies and confers protection (Cox et al.,
1977). This yields an unambiguous single-antigen ground truth (G = protective antigen) against
which a tool can be tested for correctly top-ranking the true antigen.

### A.2 Host selection (four mammalian hosts)

The viral arm's principal purpose is to exercise VaxForge's host-specific MHC epitope capability,
which the bacterial-oriented comparators do not provide. Four mammalian hosts — human, bovine,
mouse and pig — were selected because all four are genuine natural mammalian hosts of RABV, which
infects essentially all mammals (Fooks et al., 2017); each host therefore permits a biologically
honest MHC presentation rather than an artificial host assignment. Rabies specifically enables this
four-host matrix where other candidates would not: influenza is marginal in cattle, foot-and-mouth
disease virus does not infect humans, and SARS-CoV-2 is not a natural bovine/porcine host. The four
hosts are also practically meaningful — human (public-health target), cattle and pigs (veterinary
and economic burden, with pigs epidemiologically close to humans), and mouse (the standard rabies
challenge/laboratory model). The global human and livestock burden of rabies further motivates the
human and agricultural hosts (Hampson et al., 2015). We report host-specific epitope repertoires as
a unique-capability demonstration, not as a head-to-head comparison, since the comparators cannot
select a host.

### A.3 Definition of ground-truth protective antigens

For each pathogen, the ground-truth set comprises antigens with literature-documented protective or
vaccine-candidate status; every antigen and its supporting citation are listed in
`PROTECTIVE_ANTIGENS.md` (all DOIs verified 2026-09-20). Evidence strength is reported per antigen
and is not overstated: it ranges from direct active-immunization protection (e.g. *S. aureus* IsdB,
Kuklin et al., 2006; *L. monocytogenes* LLO, Phelps et al., 2020; *Salmonella* OmpD, Gil-Cruz et
al., 2009), through passive-immunization protection (e.g. *L. monocytogenes* ActA, Asano et al.,
2016), to well-established invasion virulence factors that are investigated candidate antigens
(*L. monocytogenes* InlA/InlB, Gaillard et al., 1991; Dramsi et al., 1995). This transparency is
consistent with our treatment of ActA as a legitimate hard antigen (it also ranks lowest in the
independent comparator; `BENCHMARK_SUMMARY.md` §1b).

### A.4 Competitor selection

Both comparators are justified and were run under default parameters with identical inputs (a fair
condition, no threshold tuning). NERVE 2.0 is a current, AI-assisted, locally installable reverse-
vaccinology environment (Conte et al., 2024); Vaxign-ML is an independent supervised-ML predictor of
bacterial protective antigens (Ong et al., 2020). That the independent-ML tool (Vaxign-ML) also
attains high recall shows the soft-filter behaviour distinguishing VaxForge is not a single-tool
artefact.

---

## B. Pipeline scale-safety validation

During whole-proteome runs of the extended bacterial arm we identified, isolated and corrected two
scale-dependent defects that had artefactually depressed VaxForge recall; we report them explicitly
as part of pipeline validation. Both were "single oversized batch" failures that manifested only at
full-proteome scale (initial runs: *S. aureus* 3/5, *L. monocytogenes* 3/5, *S.* Typhimurium 0/5),
while a small smoke-test subset returned 5/5 for all three species — an inconsistency that triggered
a systematic root-cause search.

**(i) Subcellular localization (PSORTb).** The localization step originally submitted every protein
to a single container call under a fixed 1800 s subprocess timeout. For the largest proteome
(4,533 proteins) the call exceeded the timeout and returned empty, causing the funnel to fall back
to heuristic localization; genuine surface antigens (e.g. outer-membrane porins) were then
mis-assigned as cytoplasmic and eliminated. The fix partitions the input into sub-batches of ≤800
proteins so each container call completes within the timeout, and merges the results.

**(ii) MHC binding (netMHCpan).** The MHC-I step originally called netMHCpan once with the peptides
of all proteins concatenated. At whole-proteome scale (~1.37 million unique 9-mers) the tool failed,
yielding near-zero coverage so that most proteins appeared to have no MHC epitopes. Because peptides
were processed in alphabetical order, which antigens survived was effectively random — accounting
for the seemingly arbitrary "losses" across species. The fix chunks peptides into blocks of 20,000
and runs the chunks in parallel across available cores (ThreadPoolExecutor); validation on a
50,528-peptide set recovered 100% coverage (e.g. OmpA 342/342 9-mers).

Control experiments established that the toxicity and allergenicity filters were **not** responsible:
their per-peptide behaviour is deterministic and identical between single-peptide and full-batch
evaluation. After both fixes, all three species were re-run under default parameters, giving the
final valid results (*S. aureus* 5/5, AUC 1.00; *L. monocytogenes* 4/5; *S.* Typhimurium 5/5,
AUC 0.99; total VaxForge 14/15). The single remaining miss, ActA, is a legitimate hard antigen that
also ranks lowest in the independent comparator (Vaxign-ML), not a residual artefact. These findings
document a scale-safe pipeline and are reported to preclude the misleading conclusion that the
initial low recall reflected the method rather than an engineering-scale defect.

---

## References (Methods; DOI-verified 2026-09-20 — see `PROTECTIVE_ANTIGENS.md` for the antigen set)

- Pizza M, et al. *Science*. 2000;287(5459):1816–1820. doi:10.1126/science.287.5459.1816
- Giuliani MM, et al. *PNAS*. 2006;103(29):10834–10839. doi:10.1073/pnas.0603940103
- Cox JH, Dietzschold B, Schneider LG. *Infect Immun*. 1977;16(3):754–759. doi:10.1128/iai.16.3.754-759.1977
- Fooks AR, et al. *Nat Rev Dis Primers*. 2017;3:17091. doi:10.1038/nrdp.2017.91
- Hampson K, et al. *PLoS Negl Trop Dis*. 2015;9(4):e0003709. doi:10.1371/journal.pntd.0003709
- Conte A, et al. *BMC Bioinformatics*. 2024;25:377. doi:10.1186/s12859-024-06004-0
- Ong E, et al. *Bioinformatics*. 2020;36(10):3185–3191. doi:10.1093/bioinformatics/btaa119
- Kuklin NA, et al. *Infect Immun*. 2006;74(4):2215–2223. doi:10.1128/iai.74.4.2215-2223.2006
- Phelps CC, et al. *Vaccine*. 2020;38(36):5803–5813. doi:10.1016/j.vaccine.2020.06.049
- Gil-Cruz C, et al. *PNAS*. 2009;106(24):9803–9808. doi:10.1073/pnas.0812431106
- Asano K, et al. *Sci Rep*. 2016;6:39628. doi:10.1038/srep39628
- Gaillard JL, et al. *Cell*. 1991;65(7):1127–1141. doi:10.1016/0092-8674(91)90009-n
- Dramsi S, et al. *Mol Microbiol*. 1995;16(2):251–261. doi:10.1111/j.1365-2958.1995.tb02297.x

> Bakteriyel ground-truth antijenlerin tam atıf listesi (18 kaynak) `PROTECTIVE_ANTIGENS.md`'de.
</content>
