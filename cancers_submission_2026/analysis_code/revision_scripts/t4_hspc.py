"""T4 (R1-7): haematopoietic stem / progenitor programs across immune ecotypes.
Modules: Azimuth 2023 human bone-marrow reference (Hao et al. 2021), retrieved via Enrichr.
Scoring matched to the existing pipeline: gseapy.ssgsea, sample_norm_method='rank', log2(TPM+1)."""
import numpy as np, pandas as pd, gseapy as gp, json
from scipy import stats

UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
WANT={ # Azimuth term -> module label
 "Bone Marrow-L2-Hematopoeitic Stem Cell":"HSC",
 "Bone Marrow-L2-Lymphoid Primed Multipotent Progenitor":"LMPP",
 "Bone Marrow-L2-Granulocyte Monocyte Progenitor":"GMP",
 "Bone Marrow-L2-Common Lymphoid Progenitor":"CLP",
 "Bone Marrow-L2-Erythroid Megakaryocyte Progenitor":"EMP (lineage control)",
 "Bone Marrow-L2-CD14 Monocyte":"CD14 monocyte (positive control)",
 "Bone Marrow-L2-Macrophage":"Macrophage (positive control)",
}
sets={}
for line in open("Azimuth_2023.gmt"):
    f=line.rstrip("\n").split("\t")
    if f[0] in WANT:
        sets[WANT[f[0]]]=[g.split(",")[0] for g in f[2:] if g.strip()]
for k,v in sets.items(): print(f"  {k:36} {len(v)} genes")

tpm=pd.read_csv(f"{UP}/tpm_for_cibersortx.tsv",sep="\t",index_col=0)
eco=pd.read_csv(f"{UP}/ecotype_LM22_main_k3_annotated.tsv",sep="\t").set_index("Kids_First_Biospecimen_ID")
tpm=tpm[[c for c in eco.index if c in tpm.columns]]
print("TPM restricted to main cohort:",tpm.shape)
for k,v in sets.items(): print(f"  {k:36} {sum(g in tpm.index for g in v)}/{len(v)} present")

log_tpm=np.log2(tpm.astype(float)+1.0)
res=gp.ssgsea(data=log_tpm,gene_sets=sets,sample_norm_method="rank",no_plot=True,
              threads=2,min_size=3,max_size=500,permutation_num=0,outdir=None)
s=res.res2d.copy(); s["NES"]=pd.to_numeric(s["NES"],errors="coerce")
s=s.pivot(index="Name",columns="Term",values="NES").astype(float)
s.index.name="Kids_First_Biospecimen_ID"
s=s.join(eco["ecotype"])
s.to_csv("T4_HSPC_ssGSEA_scores.tsv",sep="\t")
print("\nn scored:",len(s),"| ecotype counts:",s.ecotype.value_counts().to_dict())

ORDER=["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]
mods=[m for m in sets]
rows=[]
for m in mods:
    grp=[s.loc[s.ecotype==e,m].values for e in ORDER]
    H,p=stats.kruskal(*grp)
    eps=(H-len(ORDER)+1)/(len(s)-len(ORDER))
    r=dict(module=m,KW_H=round(H,2),p=p,epsilon_sq=round(eps,3))
    for e,g in zip(ORDER,grp): r[f"mean_{e}"]=round(g.mean(),4)
    rows.append(r)
kw=pd.DataFrame(rows)
kw["q_BH"]=stats.false_discovery_control(kw.p)
kw=kw.sort_values("p")
print("\n=== Kruskal-Wallis across ecotypes ===")
print(kw[["module","KW_H","p","q_BH","epsilon_sq","mean_Lymphocyte-inflamed","mean_Myeloid-dominant","mean_Immune-desert"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
kw.to_csv("T4_KW_by_ecotype.tsv",sep="\t",index=False)

# Dunn post hoc (BH within module)
def dunn(vals,labels):
    allv=np.concatenate(vals); rk=stats.rankdata(allv); ns=[len(v) for v in vals]
    idx=np.cumsum([0]+ns); N=len(allv)
    mr=[rk[idx[i]:idx[i+1]].mean() for i in range(len(vals))]
    _,cnt=np.unique(allv,return_counts=True); tie=(cnt**3-cnt).sum()
    sig2=(N*(N+1)/12)-tie/(12*(N-1))
    out=[]
    for i in range(len(vals)):
        for j in range(i+1,len(vals)):
            z=(mr[i]-mr[j])/np.sqrt(sig2*(1/ns[i]+1/ns[j]))
            out.append((labels[i],labels[j],z,2*stats.norm.sf(abs(z))))
    return out
drows=[]
for m in mods:
    vals=[s.loc[s.ecotype==e,m].values for e in ORDER]
    for a,b,z,p in dunn(vals,ORDER): drows.append(dict(module=m,group1=a,group2=b,z=round(z,2),p=p))
dn=pd.DataFrame(drows); dn["q_BH"]=stats.false_discovery_control(dn.p)
dn.to_csv("T4_Dunn_posthoc.tsv",sep="\t",index=False)
print("\n=== Dunn post hoc (BH across all 21 comparisons) ===")
print(dn.to_string(index=False,float_format=lambda x:f"{x:.3g}"))
