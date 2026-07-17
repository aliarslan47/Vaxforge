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
    ("Lokalizasyon — bakteri", "PSORTb 3.0",
     "Yu NY, Wagner JR, Laird MR, Melli G, Rey S, Lo R, et al. PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. Bioinformatics. 2010;26(13):1608-1615.",
     "10.1093/bioinformatics/btq249"),
    ("Lokalizasyon — yedek/ökaryot", "DeepLoc-2.1",
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
    # ── MEV konstruktu (fizikokimyasal + inşa) ──────────────────────────────
    ("MEV — fizikokimyasal (MW, pI, GRAVY)", "ExPASy ProtParam",
     "Gasteiger E, Hoogland C, Gattiker A, Duvaud S, Wilkins MR, Appel RD, Bairoch A. Protein Identification and Analysis Tools on the ExPASy Server. In: Walker JM, ed. The Proteomics Protocols Handbook. Humana Press; 2005:571-607. (Biopython Bio.SeqUtils.ProtParam ile hesaplandı).",
     "10.1385/1-59259-890-0:571"),
    ("MEV — hidropati (GRAVY)", "Kyte-Doolittle",
     "Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982;157(1):105-132.",
     "10.1016/0022-2836(82)90515-0"),
    ("MEV — kararlılık (instability index)", "Guruprasad et al. 1990",
     "Guruprasad K, Reddy BVB, Pandit MW. Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. Protein Eng. 1990;4(2):155-161. (indeks < 40 -> kararlı).",
     "10.1093/protein/4.2.155"),
    ("MEV — termostabilite (aliphatic index)", "Ikai 1980",
     "Ikai A. Thermostability and aliphatic index of globular proteins. J Biochem. 1980;88(6):1895-1898.",
     "10.1093/oxfordjournals.jbchem.a133168"),
    ("MEV — rijit linker (adjuvan ayırıcı)", "EAAAK linker (Arai et al. 2001)",
     "Arai R, Ueda H, Kitayama A, Kamiya N, Nagamune T. Design of the linkers which effectively separate domains of a bifunctional fusion protein. Protein Eng. 2001;14(8):529-532.",
     "10.1093/protein/14.8.529"),
    ("MEV — HTL linker (GPGPG)", "Livingston et al. 2002",
     "Livingston B, Crimi C, Newman M, Higashimoto Y, Appella E, Sidney J, Sette A. A rational strategy to design multiepitope immunogens based on multiple Th lymphocyte epitopes. J Immunol. 2002;168(11):5499-5506. (GPGPG: bağlantı immünojenitesini azaltır, HTL işleme).",
     "10.4049/jimmunol.168.11.5499"),
    ("MEV — CTL linker (AAY) ve string-of-beads", "Schubert & Kohlbacher 2016",
     "Schubert B, Kohlbacher O. Designing string-of-beads vaccines with optimal spacers. Genome Med. 2016;8:9. (AAY: proteazomal kesim; KK: B-hücre epitop salımı).",
     "10.1186/s13073-016-0263-6"),
    ("MEV — adjuvan (β-defensin, TLR agonisti)", "Biragyn et al. 2002",
     "Biragyn A, Ruffini PA, Leifer CA, Klyushnenkova E, Shakhov A, Chertov O, et al. Toll-like receptor 4-dependent activation of dendritic cells by β-defensin 2. Science. 2002;298(5595):1025-1029. (N-ucu β-defensin adjuvanı; alternatif: 50S ribozomal L7/L12).",
     "10.1126/science.1075565"),
    ("MEV — içsel düzensizlik (disorder)", "metapredict V3",
     "Emenecker RJ, Griffith D, Holehouse AS. Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure. Biophys J. 2021;120(20):4312-4319. (%disordered, aracın kendi IDR çağrılarından).",
     "10.1016/j.bpj.2021.08.039"),
    ("MEV — ikincil yapı (%helix/strand/coil)", "S4PRED (Moffat & Jones 2021)",
     "Moffat L, Jones DT. Increasing the accuracy of single sequence prediction methods using a deep semi-supervised learning framework. Bioinformatics. 2021;37(21):3744-3751. (tek-dizi C/H/E tahmini).",
     "10.1093/bioinformatics/btab491"),
    ("MEV — çözünürlük (scaled solubility)", "Protein-Sol (Hebditch et al. 2017)",
     "Hebditch M, Carballo-Amador MA, Charonis S, Curtis R, Warwicker J. Protein-Sol: a web tool for predicting protein solubility from sequence. Bioinformatics. 2017;33(19):3098-3100. (ölçekli çözünürlük; >0.45 çözünür).",
     "10.1093/bioinformatics/btx345"),
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
