# Vaxforge

Web arayüzlü, ajan-destekli, in silico **ters aşılama (reverse vaccinology)** pipeline'ı — bir patojen dizisinden sıralanmış çok-epitoplu mRNA aşı yapısına ve tam atıflı bir rapora.

[![type](https://img.shields.io/badge/type-reverse%20vaccinology-0d6b8f)](https://github.com/aliarslan47/Vaxforge)
[![interface](https://img.shields.io/badge/interface-Streamlit%20web-2f8f5b)](https://github.com/aliarslan47/Vaxforge)
[![tools](https://img.shields.io/badge/tools-real%20%C2%B7%20cited-c07211)](https://github.com/aliarslan47/Vaxforge)

**Türkçe** · [English](README.md)

## Nedir?

Vaxforge, Forge ailesinin aşı-tasarımı üyesidir — bir patojen dosyasını sıralanmış aşı adaylarına çeviren bir ters aşılama pipeline'ı. Biyolojiyi deterministik bir bilimsel çekirdek yapar; LLM sarmalayıcı yalnızca planlar, yorumlar ve raporlar. Bir Streamlit web arayüzüyle gelir.

## Ne yapar?

Bir patojen dosyası (FASTA/FASTQ) yükle → sistem girdiyi otomatik algılar → virülans faktörlerini ve aşı hedeflerini tarar → antijen hunisini çalıştırır → seçilebilir konakçılarda B- ve T-hücre epitoplarını tahmin eder → hayatta kalma filtreleri (toksisite, alerjenite) → adayları yapılandırılabilir bir adaylık skoruyla sıralar → çok-epitoplu bir mRNA yapısı kurar → tam araç atıflı, yayın-tarzı bir rapor üretir.

- **Gerçek, yayınlanmış araçlar**; her biri raporda (PDF/HTML/JSON) ve arayüzde atıflı.
- **Konakçı-seçilebilir, çok-organizmalı MHC** (insan, fare, sığır, domuz, tavuk) → konakçı-başına sunum haritası.
- **Sihirli sayı yok**: eşikler `config/thresholds.yaml`'da (organizma ön-ayarları), arayüzden düzenlenebilir ve her rapora yazılır.
- **Dürüst etiketleme**: gerçek bir araç yoksa açıkça-etiketli bir fallback kullanılır.

## Kurulum

```bash
git clone https://github.com/aliarslan47/Vaxforge.git
cd Vaxforge
python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt
```

Harici araçlar `tools/` altındadır (git-ignore). Ücretsiz araçlar otomatik kurulur (DIAMOND, VFDB, Swiss-Prot, ToxinPred2, IApred, DeepLoc-2.1); lisanslı DTU HealthTech araçları (NetMHCpan, NetMHCIIpan, SignalP-5.0, TMHMM-2.0, BepiPred-1.0) manuel akademik indirme gerektirir — bkz. `tools/README.md`. Bir araç eksikse pipeline etiketli bir fallback ile yine çalışır.

## Kullanım

```bash
streamlit run app.py
```

Bir FASTA/FASTQ dosyası sürükle (ya da örnek seç), patojen profilini ve konakçı(lar)ı seç, eşikleri gözden geçir ve çalıştır. Çıktılar: yayın-tarzı **PDF**, **HTML** panosu, sıralı **CSV**, peptit **FASTA**, **GenBank** mRNA yapısı ve tam-run **JSON**.

## Modüller

Bilimsel çekirdek gerçek, atıflı araçlarla çalışır; GPU-bağımlı yapısal adımlar ertelenmiştir.

| Adım | Araç / yöntem | Durum |
|---|---|---|
| Girdi otomatik algılama | FASTA/FASTQ · nt/protein · genom/CDS/okuma | yerleşik |
| Keşif (virülans faktörleri) | DIAMOND + VFDB | ✅ gerçek |
| Antijen hunisi — lokalizasyon | DeepLoc-2.1 | ✅ gerçek |
| Antijen hunisi — transmembran | TMHMM-2.0 | ✅ gerçek |
| Antijen hunisi — sinyal peptidi | SignalP-5.0 | ✅ gerçek |
| Antijen hunisi — antijenite | IApred | ✅ gerçek |
| Antijen hunisi — konakçı homolojisi (güvenlik) | DIAMOND vs insan Swiss-Prot | ✅ gerçek |
| Epitoplar — B-hücre | BepiPred-1.0 | ✅ gerçek |
| Epitoplar — MHC-I / MHC-II | NetMHCpan / NetMHCIIpan (yerel + IEDB) | ✅ gerçek |
| Hayatta kalma — toksisite / alerjenite | ToxinPred2 · FAO/WHO 6-mer + UniProt | ✅ gerçek |
| Adaylık skorlaması | ağırlıklı, yapılandırılabilir | ✅ gerçek |
| mRNA yapısı | linker + adjuvan + insan kodon-opt + GC/CAI | ✅ gerçek |
| Peptit–MHC yapısı + MD | AlphaFold + docking/MD | ⏸️ ertelendi (GPU) |

Tam atıflar (`vaxforge/citations.py`) ve proje yapısı depoda.

---

Forge ailesi: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bakteri) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virüs/faj) · [MicrobiomeForge](https://github.com/aliarslan47/MicrobiomeForge) (mikrobiyom) · **Vaxforge** (ters aşılama) · [ImmForge](https://github.com/aliarslan47/ImmForge) (bağışıklık simülasyonu) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG üreticisi).
