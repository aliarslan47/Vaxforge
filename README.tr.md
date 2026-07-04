# VaxForge 🧬

[🌐 English](README.md) | 🇹🇷 Türkçe

Ajan destekli **in silico reverse vaccinology** pipeline'ı ve arayüzü (prototip).

Bir patojen dosyası (FASTA/FASTQ) yükle → sistem tipini **otomatik tanısın** →
ne yapacağını **planlasın** → virülans/aşı hedeflerini, B/T-hücre epitoplarını
çıkarsın → en güçlü peptit adaylarını **adaylık puanıyla** sıralasın → çok-epitoplu
**mRNA konstrüktü** kursun → rapor + veri paketi üretsin.

## Durum (bu sürüm)

Hafif omurganın ilk dilimi hazır ve çalışıyor:

- ✅ **Dosya tanıma** — FASTA/FASTQ, nükleotid/protein, molekül türü (reads/genome/cds/proteome), `.gz` desteği
- ✅ **Pipeline planlama** — girdiye göre adım listesi; GPU'lu ağır adımlar (AlphaFold peptit-MHC, docking/MD) `deferred` işaretli
- ✅ **Merkezi eşik/config** — her aracın cut-off'u `config/thresholds.yaml`'de; organizmaya özel presetler (bakteri/virüs/parazit); kodda magic number yok
- ✅ **Streamlit arayüzü** — yükle → tanı → plan → düzenlenebilir eşikler
- ⏳ Sıradaki: keşif (küratörlü DB), antijen hunisi, epitop tahmini, sağ kalım elemesi, skorlama, mRNA, rapor modülleri

## Kurulum & çalıştırma

```bash
cd /home/msi-nb/vaxforge
python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Tarayıcıda açılan arayüze `data/samples/` içindeki örnek dosyalardan birini
yükleyerek tanıma + plan + eşik ekranını görebilirsin.

## Mimari kararlar

- **Deterministik hat + LLM sarmalayıcı:** Bilimsel adımlar sabit/tekrarlanabilir; LLM yalnızca yorum ve raporlama.
- **İki katmanlı çalışma:** Bu makine (GPU'suz) = geliştirme + hafif omurga; ağır modüller ayrı GPU makinesinde/bulutta.
- **Konaktan bağımsız ama çok-organizmalı MHC:** Peptitler birden çok türün MHC'sine karşı test edilir (tür-kapsamı haritası).
- **Eşikler koda gömülmez:** Hepsi config'te, organizmaya göre, arayüzden düzenlenebilir, rapora yazılır.

## Yapı

```
config/thresholds.yaml     # tüm eşikler + organizma presetleri
vaxforge/detect.py         # dosya tanıma
vaxforge/plan.py           # tanımaya göre pipeline planı
vaxforge/config_loader.py  # eşik yükleme/doğrulama/çözme
app.py                     # Streamlit arayüzü
data/samples/              # örnek FASTA/FASTQ
```
