"""Yeni 6 bakterinin UniProt referans proteomunu indir + genom boyutunu NCBI'dan al.
Hepsi doğrulanmış aksesyon (uydurma yok). Çıktı: data/validation/<slug>_proteome.faa
"""
import urllib.request, urllib.parse, json, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "validation"

# slug, upid, organism, taxid_strain, taxid_species, assembly
BUGS = [
    ("coxiella",   "UP000002671", "Coxiella burnetii RSA 493 / Nine Mile phase I", 227377, 777,   "GCA_000007765.2", "negative"),
    ("ssuis",      "UP000002055", "Streptococcus suis P1/7",                       218494, 1307,  "GCA_000091905.1", "positive"),
    ("brucella",   "UP000007104", "Brucella suis biovar 1 strain 1330",            204722, 29461, "GCA_000223195.1", "negative"),
    ("mbovis",     "UP000001419", "Mycobacterium bovis AF2122/97",                 233413, 1765,  "GCA_000195835.3", "positive"),
    ("leptospira", "UP000001408", "Leptospira interrogans serovar Lai 56601",      189518, 173,   "GCA_000092565.1", "negative"),
    ("banthracis", "UP000000427", "Bacillus anthracis Ames Ancestor",              261594, 1392,  "GCA_000007845.1", "positive"),
]

def genome_mb(acc):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/dataset_report"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
        rep = d["reports"][0]
        tl = int(rep["assembly_stats"]["total_sequence_length"])
        return round(tl/1e6, 2)
    except Exception as e:
        return f"NA({e})"

def download_proteome(upid, out):
    url = (f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{upid}"
           "&format=fasta&includeIsoform=false")
    with urllib.request.urlopen(url, timeout=300) as r:
        data = r.read()
    out.write_bytes(data)
    n = data.count(b">")
    res = sum(len(l) for l in data.split(b"\n") if l and not l.startswith(b">"))
    return n, res

rows = []
for slug, upid, org, ts, tsp, asm, gram in BUGS:
    mb = genome_mb(asm)
    out = OUT / f"{slug}_proteome.faa"
    n, res = download_proteome(upid, out)
    rows.append((mb if isinstance(mb,float) else 9999, slug, org, upid, n, res, mb, gram, ts, tsp, asm))
    print(f"{slug:11s} {upid} | {n:5d} protein | {res:8d} res | genom {mb} Mb | {org}")
    time.sleep(0.5)

print("\n===== GENOM KÜÇÜKTEN BÜYÜĞE SIRA =====")
for i,(k,slug,org,upid,n,res,mb,gram,ts,tsp,asm) in enumerate(sorted(rows),1):
    print(f"  {i}. {slug:11s} {mb} Mb · {n} protein · {gram} · {org}")
