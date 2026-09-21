"""Stream only the genes we need out of the two large TPM matrices (genes = rows)."""
import gzip, csv, sys, numpy as np, pandas as pd
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"

# signature genes
sigs={}
for line in open(f"{UP}/brain_immune_signatures.gmt"):
    f=line.rstrip("\n").split("\t"); sigs[f[0]]=[g for g in f[2:] if g.strip()]
sig_genes=set(g for v in sigs.values() for g in v)
# ecotype-defining genes: top 300 by q from each adjusted contrast
deg=pd.read_csv(f"{UP}/adjusted_DEG_all_contrasts.tsv",sep="\t")
top=set()
for c,g in deg.groupby("contrast"):
    up=g[(g.q_BH<0.05)&(g.adjusted_log2FC>1)].nsmallest(300,"q_BH")["gene"]
    top|=set(up)
# canonical identity markers for verification
MARK=["PTPRC","CD3D","CD3E","CD2","CD8A","CD4","IL7R","TRAC","CSF1R","AIF1","C1QA","C1QB","CD14",
      "TYROBP","P2RY12","CX3CR1","ITGAM","CD68","MBP","PLP1","SOX2","OLIG1","OLIG2","PDGFRA","GFAP","EGFR"]
WANT=sig_genes|top|set(MARK)
print(f"signature genes {len(sig_genes)} | ecotype DEG top {len(top)} | total wanted {len(WANT)}")
json_out={"n_sig":len(sig_genes),"n_deg":len(top),"n_want":len(WANT)}

def extract(path,out):
    with gzip.open(path,"rt",newline="") as f:
        r=csv.reader(f)
        hdr=next(r); cells=hdr[1:]
        rows=[]; names=[]
        for i,row in enumerate(r):
            if row[0] in WANT:
                names.append(row[0]); rows.append(np.asarray(row[1:],dtype=np.float32))
            if i%5000==0: print(f"  {path.split('/')[-1]}: {i} rows, kept {len(names)}",flush=True)
        M=pd.DataFrame(np.vstack(rows),index=names,columns=cells)
        M=M[~M.index.duplicated()]
        M.to_pickle(out); print(f"  -> {out} {M.shape}")
        return M

extract("data/immune_tenx_tpm.csv.gz","w2_immune_tenx.pkl")
extract("data/tumor_ss2_tpm.csv.gz","w2_tumor_ss2.pkl")
import json; json.dump(json_out,open("w2_genesets_meta.json","w"))
pd.Series({k:";".join(v) for k,v in sigs.items()}).to_pickle("w2_sigs.pkl")
pd.Series(sorted(top)).to_pickle("w2_topdeg.pkl")
