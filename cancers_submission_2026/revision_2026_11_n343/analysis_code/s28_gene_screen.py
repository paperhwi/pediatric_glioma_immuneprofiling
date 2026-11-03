"""Per-gene Cox screen of the shared immune-presence genes (Figure 6D): OS ~ gene (z of log2 TPM+1 across n = 343) +
molecular group + age + sex; BH across genes (as Rev344_06 supp.py S31)."""
import sys, glob; sys.path.insert(0,"/home/claude/q343/src"); from common import *
a=load(); s=a[a.os_eval]; tpm=pd.read_pickle(ROOT+"tpm343.pkl")
sh=pd.read_csv(glob.glob(ROOT+"out/T_shared_genes_*.tsv")[0],sep="\t").gene.tolist()
x=np.log2(tpm.loc[[g for g in sh if g in tpm.index]]+1); z=x.sub(x.mean(1),axis=0).div(x.std(1)+1e-9,axis=0)
out=[]
for g in z.index:
    c=cox(design(s,extra=z.loc[[g]].T.rename(columns={g:"gene"}),eco=False)); r=c.summary.loc["gene"]
    out.append(dict(gene=g,HR_per_SD=r["exp(coef)"],lo=r["exp(coef) lower 95%"],hi=r["exp(coef) upper 95%"],p=r["p"]))
T=pd.DataFrame(out); T["q_BH"]=bh(T.p); T=T.sort_values("p"); T.to_csv(ROOT+"out/T_gene_screen_shared_genes.tsv",sep="\t",index=False)
print(len(T),"q<0.05:",int((T.q_BH<0.05).sum()),"protective:",int(((T.q_BH<0.05)&(T.HR_per_SD<1)).sum()))
