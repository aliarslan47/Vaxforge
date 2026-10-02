"""Her bakteri için koruyucu antijen UniProt aksesyonlarını proteom İÇİNDEN çöz (uydurma yok).
gene araması -> protein-adı araması fallback. Çıktı: aksesyon + gen + ad + proteomda-var mı.
"""
import urllib.request, urllib.parse, json, time
# slug, upid, {etiket: [gen adayları / protein-adı anahtarları]}
ANTIGENS = {
 "coxiella":  ("UP000002671", {
    "Com1":"com1", "Mip":"mip", "OmpH":"ompH", "GroEL":"groEL"}),
 "ssuis":     ("UP000009077", {
    "Suilysin(SLY)":"sly", "MRP":"mrp", "Enolase":"eno", "GAPDH":"gap", "IgdE/ Meningitis":"srtA"}),
 "brucella":  ("UP000007104", {
    "Omp25":"omp25", "Omp31":"omp31", "BP26":"bp26", "SodC":"sodC", "L7/L12(RplL)":"rplL"}),
 "mbovis":    ("UP000001419", {
    "ESAT-6(EsxA)":"esxA", "CFP-10(EsxB)":"esxB", "Ag85B(FbpB)":"fbpB", "MPT64":"mpt64", "HspX(Acr)":"hspX"}),
 "leptospira":("UP000001408", {
    "LigA":"ligA", "LigB":"ligB", "LipL32":"lipL32", "LipL41":"lipL41", "OmpL1":"ompL1"}),
 "banthracis":("UP000000427", {
    "PA(pagA)":"pagA", "LF(lef)":"lef", "EF(cya)":"cya", "Sterne-BslA":"bslA"}),
}
def search(query):
    url=("https://rest.uniprot.org/uniprotkb/search?query="+urllib.parse.quote(query)
         +"&format=tsv&fields=accession,gene_names,protein_name,length&size=5")
    with urllib.request.urlopen(url,timeout=60) as r: return r.read().decode().strip().splitlines()
for slug,(upid,ags) in ANTIGENS.items():
    print(f"\n===== {slug} ({upid}) =====")
    for label,gene in ags.items():
        rows=search(f"proteome:{upid} AND gene:{gene}")
        if len(rows)<=1:
            rows=search(f'proteome:{upid} AND protein_name:"{label.split("(")[0]}"')
        if len(rows)<=1:
            print(f"  {label:18s} -> BULUNAMADI (gene:{gene})")
        else:
            acc,gn,pn,ln=rows[1].split("\t")+["","",""] if len(rows[1].split("\t"))<4 else rows[1].split("\t")[:4]
            print(f"  {label:18s} -> {acc:10s} | {gn:14s} | {pn[:45]} | {ln}aa")
        time.sleep(0.3)
