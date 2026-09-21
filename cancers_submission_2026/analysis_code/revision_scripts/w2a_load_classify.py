"""W2.1 — load GSE102130 (Filbin 2018 H3K27M glioma, Smart-seq2) and classify cells by markers."""
import numpy as np, pandas as pd
D=pd.read_csv("data/GSE102130.txt.gz",sep="\t",index_col=0)
D.index=D.index.astype(str)
print("raw matrix:",D.shape)
cols=pd.Series(D.columns)
pat=cols.str.split("-").str[0]
DERIV=["BCH869_DGC100","BCH869_DGC75","BCH869_GF","BCH869_GS","BCH869_PDX","BCH869_SF"]
keep=~pat.isin(DERIV+["Oligo"])
print("primary-tumour cells:",int(keep.sum()),"| excluded culture/PDX:",int(pat.isin(DERIV).sum()),
      "| normal oligo reference:",int((pat=='Oligo').sum()))
E=D.loc[:,keep.values].astype(np.float32)
patient=pat[keep.values].values
# values are already log2(TPM/10+1) in Filbin's processed table? check scale
print("value range:",float(E.values.min()),float(np.percentile(E.values,99.9)),float(E.values.max()))
