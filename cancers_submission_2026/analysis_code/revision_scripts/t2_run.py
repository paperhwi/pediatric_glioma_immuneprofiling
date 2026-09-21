import numpy as np, pandas as pd
from scipy import stats
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
d=pd.read_csv(f"{UP}/adjusted_DEG_all_contrasts.tsv",sep="\t")
A="Lymphocyte_inflamed_vs_Immune_desert"; B="Myeloid_dominant_vs_Immune_desert"; C="Lymphocyte_inflamed_vs_Myeloid_dominant"
w=d.pivot(index="gene",columns="contrast",values="adjusted_log2FC")
q=d.pivot(index="gene",columns="contrast",values="q_BH")
sig=lambda c,dirn=+1:set(w.index[(q[c]<0.05)&((w[c]>1) if dirn>0 else (w[c]<-1))])
up={k:sig(k) for k in (A,B,C)}; dn={k:sig(k,-1) for k in (A,B,C)}
print("=== UP (q<0.05, adj log2FC > 1) ===")
for k in (A,B,C): print(f"  {k:42} {len(up[k]):5}   DOWN {len(dn[k])}")

sA,sB=up[A],up[B]
print(f"\nShared by both vs-Immune-desert contrasts : {len(sA&sB):5}")
print(f"Unique to Lymphocyte-inflamed vs Desert   : {len(sA-sB):5}")
print(f"Unique to Myeloid-dominant vs Desert      : {len(sB-sA):5}")
print(f"Jaccard(LI vs ID, MD vs ID)               : {len(sA&sB)/len(sA|sB):.3f}")
print(f"Of the {len(sB)} MD-vs-Desert genes, {len(sA&sB)/len(sB)*100:.1f}% are also LI-vs-Desert")

# is the shared set explained by the LI>MD gradient?
sh=sorted(sA&sB)
sub=w.loc[sh]
print(f"\nAmong the {len(sh)} shared genes: mean adj log2FC  LI vs ID = {sub[A].mean():.3f},  MD vs ID = {sub[B].mean():.3f}")
print(f"  paired Wilcoxon LI-vs-ID > MD-vs-ID : p = {stats.wilcoxon(sub[A],sub[B],alternative='greater').pvalue:.3e}")
print(f"  median ratio of effect sizes        : {np.median(sub[A]/sub[B]):.2f}x")
# orthogonal (Deming) regression slope on ALL genes significant in either
both=sorted(sA|sB); x=w.loc[both,B].values; y=w.loc[both,A].values
xc,yc=x-x.mean(),y-y.mean()
sxx,syy,sxy=(xc**2).mean(),(yc**2).mean(),(xc*yc).mean()
slope=((syy-sxx)+np.sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy)
r=np.corrcoef(x,y)[0,1]
print(f"\nOrthogonal-regression slope (LI-vs-ID on MD-vs-ID), n={len(both)}: {slope:.2f}   Pearson r = {r:.3f}")

# the discriminating contrast
print(f"\nLI vs MD (the contrast that actually separates the two non-desert ecotypes): "
      f"{len(up[C])} up / {len(dn[C])} down")
topC=d[(d.contrast==C)&(d.q_BH<0.05)&(d.adjusted_log2FC>1)].nsmallest(15,"q_BH")["gene"].tolist()
print("  top 15 by q:",", ".join(topC))
print(f"  of these, {sum(g in (sA&sB) for g in topC)}/15 are in the shared vs-Desert set")

rows=[dict(set_name="Shared: LI-vs-Desert AND MD-vs-Desert",n=len(sA&sB)),
      dict(set_name="Unique: LI-vs-Desert only",n=len(sA-sB)),
      dict(set_name="Unique: MD-vs-Desert only",n=len(sB-sA)),
      dict(set_name="LI-vs-Myeloid (discriminating contrast)",n=len(up[C]))]
pd.DataFrame(rows).to_csv("T2_contrast_overlap_summary.tsv",sep="\t",index=False)
pd.DataFrame({"gene":sh,"adj_log2FC_LI_vs_Desert":w.loc[sh,A].round(4),
              "adj_log2FC_MD_vs_Desert":w.loc[sh,B].round(4),
              "adj_log2FC_LI_vs_Myeloid":w.loc[sh,C].round(4)}).to_csv(
    "T2_shared_immune_presence_genes.tsv",sep="\t",index=False)
w.to_pickle("T2_w.pkl"); q.to_pickle("T2_q.pkl")
np.save("T2_sets.npy",np.array([len(sA&sB),len(sA-sB),len(sB-sA)]))
