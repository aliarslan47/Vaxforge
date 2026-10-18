<!-- Not: tools/ .gitignore kapsamında olduğu için bu kayıt benchmark/ altında tutulur. -->
# Kaldırılan araç dosyaları — 2026-10-18

Çalışma zamanında KULLANILMAYAN dosyalar silindi; **hiçbir aracın kendisi
kaldırılmadı**, yalnız eğitim verisi, yabancı-mimari derlemeleri, yinelenen
dağıtım kopyaları ve makale figürleri. Toplam kazanç **239 MB**
(tools/ 2.2 GB → 2.0 GB).

Silmeden önce iki kontrol yapıldı:
1. Depo genelinde (tools/ hariç) ve tool betiklerinin kendi içinde referans
   araması — hiçbirine atıf bulunamadı.
2. Dosyalar önce karantinaya alındı, etkilenen dört araç (IApred, netCTL,
   ToxinPred, BepiPred-3.0) gerçek diziyle koşuldu, sonra karantina boşaltıldı.

| yol | boyut | neden gereksiz |
|---|---|---|
| `IAPred/IApred-training/` | 122 MB | eğitim verisi + referans proteom; tahmin `models/IApred_SVM.joblib` ile yapılır |
| `bepipred-3.0/misc/` | 39 MB | BioLib paketinin yinelenen kopyası + eski koşum çıktıları (`output_directory`, `esm_encodings`); üretimde kök dizindeki `bepipred3_CLI.py` + `bp3/BP3Models` kullanılıyor |
| `bepipred-3.0/pip/` | 32 MB | indirme paketinin yinelenen kopyası (`BP3`, `HealtechServerDownload`) |
| `netCTL-1.2a/Linux_ia64/` | 25 MB | Itanium derlemesi; `netCTL` betiği `uname -m` ile platform seçer → bu makinede yalnız `Linux_x86_64` kullanılır |
| `netCTL-1.2a/Linux_i386/` | 4.2 MB | 32-bit derleme, aynı gerekçe |
| `netCTL-1.2a/Linux_i{486,586,686}/` | 0 | boş dizinler |
| `toxinpred2/datasets/` | 17 MB | eğitim veri setleri; tahmin `RF_model/` (228 MB, DURUYOR) ile yapılır |
| `bepipred-3.0/*.png` (3 dosya) | 2.9 MB | makale figürleri / grafiksel özet |

## Kasıtlı olarak DURAN büyük dosyalar

| yol | boyut | neden gerekli |
|---|---|---|
| `s4pred/weights/` | 428 MB | `vaxforge/secstruct.py` ikincil yapı tahmini |
| `db/` | 378 MB | insan/fare/sığır 9-mer mimikri tabloları + DIAMOND indeksleri + IEDB |
| `toxinpred2/RF_model/` | 228 MB | ToxinPred2 Model-2 hibrit tahmini |
| `netMHCpan-4.2/data`, `netMHCIIpan-4.3/data` | 190 + 125 MB | eğitilmiş allel ağları |
| `bepipred-1.0/` | 11 MB | ölü sanılmıştı; `vaxforge/immforge_score.py` kullanıyor |

## Geri getirme

Hiçbiri üretilmiş veri değil, hepsi kaynağından yeniden indirilebilir:
IAPred eğitim verisi ve ToxinPred2 datasets ilgili projelerin dağıtım
arşivlerinde; netCTL-1.2a yabancı-mimari derlemeleri DTU netCTL dağıtımında;
BepiPred-3.0 `misc/` ve `pip/` kopyaları zaten kök dizindeki kurulumun
aynısıdır (yerel olarak `bp3/` + `bepipred3_CLI.py` ile yeniden üretilebilir).
