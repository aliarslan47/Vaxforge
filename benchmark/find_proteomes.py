import urllib.request, urllib.parse, time
TARGETS = [
    ("coxiella",   "Coxiella burnetii",      "RSA 493"),
    ("ssuis",      "Streptococcus suis",     "P1/7"),
    ("brucella",   "Brucella suis",          "1330"),
    ("mbovis",     "Mycobacterium bovis",    "AF2122/97"),
    ("leptospira", "Leptospira interrogans", "Copenhageni"),
    ("banthracis", "Bacillus anthracis",     "Ames"),
]
FIELDS = "upid,organism,organism_id,protein_count,genome_assembly"
def q(name, ref_only=True):
    query = f'"{name}"' + (" AND proteome_type:1" if ref_only else "")
    url = ("https://rest.uniprot.org/proteomes/search?query="
           + urllib.parse.quote(query) + f"&format=tsv&fields={FIELDS}&size=15")
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode()
for slug, name, hint in TARGETS:
    print(f"\n===== {slug}: {name}  (REFERANS proteomlar; suş ipucu: {hint}) =====")
    try:
        txt = q(name, True).strip().splitlines()
    except Exception as e:
        print("  HATA(ref):", e); txt=[]
    if len(txt) <= 1:
        print("  (referans yok, tüm proteomlar deneniyor)")
        try: txt = q(name, False).strip().splitlines()
        except Exception as e: print("  HATA:", e); continue
    if not txt: print("  sonuç yok"); continue
    print("  ", txt[0])
    rows = txt[1:]
    rows.sort(key=lambda l: (hint.lower() not in l.lower()))
    for l in rows[:6]: print("  ", l)
    time.sleep(0.5)
