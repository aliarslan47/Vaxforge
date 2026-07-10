"""Kullanılan tüm araçların/yöntemlerin akademik atıfları.

Tek kaynak: hem PDF raporun 'Referanslar' bölümü hem pipeline özeti/JSON çıktısı
buradan beslenir. Her giriş: adım, araç, tam atıf, DOI/URL.
"""

from __future__ import annotations

# (step, tool, citation, doi_or_url)
CITATIONS = [
    ("Keşif — hizalama", "DIAMOND",
     "Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods. 2021;18:366-368.",
     "10.1038/s41592-021-01101-x"),
    ("Keşif — virülans DB", "VFDB",
     "Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. Nucleic Acids Res. 2022;50(D1):D912-D917.",
     "10.1093/nar/gkab1107"),
    ("Lokalizasyon", "DeepLoc-2.1",
     "Thumuluri V, Almagro Armenteros JJ, Johansen AR, Nielsen H, Winther O. DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. Nucleic Acids Res. 2022;50(W1):W228-W234.",
     "10.1093/nar/gkac278"),
    ("Transmembran", "TMHMM-2.0",
     "Krogh A, Larsson B, von Heijne G, Sonnhammer ELL. Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. J Mol Biol. 2001;305(3):567-580.",
     "10.1006/jmbi.2000.4315"),
    ("Sinyal peptidi", "SignalP-5.0",
     "Almagro Armenteros JJ, Tsirigos KD, Sønderby CK, et al. SignalP 5.0 improves signal peptide predictions using deep neural networks. Nat Biotechnol. 2019;37:420-423.",
     "10.1038/s41587-019-0036-z"),
    ("Antijenite", "IApred",
     "Miles S, Menafra G, Iriarte A, Chabalgoity JA. IApred: A versatile open-source tool for predicting protein antigenicity across diverse pathogens. ImmunoInformatics. 2025.",
     "10.1016/j.immuno.2025.100014"),
    ("İnsan homoloji", "UniProt/DIAMOND",
     "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res. 2023;51(D1):D523-D531. (insan Swiss-Prot; DIAMOND blastp ile).",
     "10.1093/nar/gkac1052"),
    ("MHC-I bağlanma", "NetMHCpan-4.1/4.2",
     "Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res. 2020;48(W1):W449-W454.",
     "10.1093/nar/gkaa379"),
    ("MHC-II bağlanma", "NetMHCIIpan-4.x",
     "Reynisson B, Barra C, Kaabinejadian S, Hildebrand WH, Peters B, Nielsen M. Improved prediction of MHC II antigen presentation through integration and motif deconvolution of MS MHC eluted ligand data. J Proteome Res. 2020;19(6):2304-2315.",
     "10.1021/acs.jproteome.9b00874"),
    ("IFN-γ indükleme (MHC-II)", "IFNepitope2",
     "Dhall A, Patyal S, Raghava GPS. A hybrid method for discovering interferon-gamma inducing peptides in human and mouse. Sci Rep. 2024;14(1):26859.",
     "10.1038/s41598-024-77957-8"),
    ("İmmünojenite (MHC-I)", "IEDB Class-I Immunogenicity",
     "Calis JJA, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Keşmir C, Peters B. Properties of MHC class I presented peptides that enhance immunogenicity. PLoS Comput Biol. 2013;9(10):e1003266.",
     "10.1371/journal.pcbi.1003266"),
    ("MHC anchor/cep motifi", "Anchor kalıntı kavramı + SYFPEITHI",
     "Falk K, Rötzschke O, Stevanović S, Jung G, Rammensee HG. Allele-specific motifs revealed by sequencing of self-peptides eluted from MHC molecules. Nature. 1991;351:290-296. Saper MA, Bjorkman PJ, Wiley DC. Refined structure of HLA-A2 at 2.6 Å. J Mol Biol. 1991;219:277-319. Rammensee HG, Bachmann J, Emmerich NP, Bachor OA, Stevanović S. SYFPEITHI: database for MHC ligands and peptide motifs. Immunogenetics. 1999;50:213-219. (allel motifleri NetMHCpan taramasından AMPİRİK çıkarıldı).",
     "10.1038/351290a0"),
    ("Literatür/bilinen-epitop eşleşmesi", "IEDB (Immune Epitope Database)",
     "Vita R, Blazeska N, Marrama D, et al. The Immune Epitope Database (IEDB): 2024 update. Nucleic Acids Res. 2025;53(D1):D436-D443. (aday peptitler IQ-API üzerinden deneysel doğrulanmış epitoplarla deterministik dizi eşleştirmesiyle karşılaştırıldı; LLM madenciliği kullanılmadı).",
     "10.1093/nar/gkae1092"),
    ("Popülasyon kapsamı", "IEDB Population Coverage",
     "Bui HH, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A. Predicting population coverage of T-cell epitope-based diagnostics and vaccines. BMC Bioinformatics. 2006;7:153. (allelefrequencies.net insan HLA frekansları).",
     "10.1186/1471-2105-7-153"),
    ("Antijen işleme (MHC-I)", "NetCTL-1.2 (kesim + TAP)",
     "Larsen MV, Lundegaard C, Lund O, Nielsen M. An integrative approach to CTL epitope prediction: a combined algorithm integrating MHC class I binding, TAP transport efficiency, and proteasomal cleavage predictions. Eur J Immunol. 2005;35(8):2295-2303. (yalnız allel-bağımsız kesim + TAP kullanıldı).",
     "10.1002/eji.200425811"),
    ("B-hücre epitopu", "BepiPred-3.0",
     "Clifford JN, Høie MH, Deleuran S, Peters B, Nielsen M, Marcatili P. BepiPred-3.0: Improved B-cell epitope prediction using protein language models. Protein Sci. 2022;31(12):e4497.",
     "10.1002/pro.4497"),
    ("B-hücre epitopu (yedek)", "BepiPred-1.0",
     "Larsen JEP, Lund O, Nielsen M. Improved method for predicting linear B-cell epitopes. Immunome Res. 2006;2:2.",
     "10.1186/1745-7580-2-2"),
    ("Toksisite", "ToxinPred2",
     "Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method for predicting toxicity of proteins. Brief Bioinform. 2022;23(5):bbac174.",
     "10.1093/bib/bbac174"),
    ("Alerjenite", "FAO/WHO ölçütü",
     "FAO/WHO. Evaluation of Allergenicity of Genetically Modified Foods. Report of a Joint FAO/WHO Expert Consultation. 2001. (bilinen allerjenle ≥6 ardışık aa; UniProt allergen KW-0020).",
     "https://www.fao.org/3/y0820e/y0820e00.htm"),
    ("Antijenite — yedek yöntem", "z-descriptor ACC",
     "Doytchinova IA, Flower DR. VaxiJen: a server for prediction of protective antigens... BMC Bioinformatics. 2007;8:4. (ACC özellik yöntemi; resmi model kullanılmadı). Sandberg M, et al. J Med Chem. 1998;41:2481-2491 (z-skalaları).",
     "10.1186/1471-2105-8-4"),
    ("B-hücre — yedek yöntem", "Kolaskar-Tongaonkar / Parker",
     "Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. FEBS Lett. 1990;276:172-174. Parker JMR, Guo D, Hodges RS. Biochemistry. 1986;25:5425-5432.",
     "10.1016/0014-5793(90)80535-Q"),
]


def for_report() -> list[dict]:
    return [{"step": s, "tool": t, "citation": c, "doi": d} for s, t, c, d in CITATIONS]


def as_text() -> str:
    """Pipeline özeti / konsol için numaralı düz metin."""
    lines = []
    for i, (s, t, c, d) in enumerate(CITATIONS, 1):
        ref = "https://doi.org/" + d if not d.startswith("http") else d
        lines.append(f"[{i}] {t} ({s}): {c} {ref}")
    return "\n".join(lines)
