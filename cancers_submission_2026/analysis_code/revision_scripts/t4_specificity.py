"""T4 specificity: is the GMP signal independent of mature monocyte/macrophage content?"""
import numpy as np, pandas as pd
from scipy import stats
sets={}
WANT={"Bone Marrow-L2-Hematopoeitic Stem Cell":"HSC","Bone Marrow-L2-Granulocyte Monocyte Progenitor":"GMP",
      "Bone Marrow-L2-CD14 Monocyte":"CD14 monocyte (positive control)","Bone Marrow-L2-Macrophage":"Macrophage (positive control)",
      "Bone Marrow-L2-Lymphoid Primed Multipotent Progenitor":"LMPP","Bone Marrow-L2-Common Lymphoid Progenitor":"CLP",
      "Bone Marrow-L2-Erythroid Megakaryocyte Progenitor":"EMP (lineage control)"}
for line in open("Azimuth_2023.gmt"):
    f=line.rstrip("\n").split("\t")
    if f[0] in WANT: sets[WANT[f[0]]]=[g.split(",")[0] for g in f[2:] if g.strip()]
print("=== module gene lists (Azimuth 2023 bone-marrow reference) ===")
for k in ["HSC","LMPP","GMP","CLP","EMP (lineage control)","CD14 monocyte (positive control)","Macrophage (positive control)"]:
    print(f"  {k:36} {', '.join(sets[k])}")
print("\ngene overlap GMP vs CD14 monocyte:",set(sets["GMP"])&set(sets["CD14 monocyte (positive control)"]) or "none")
print("gene overlap GMP vs Macrophage   :",set(sets["GMP"])&set(sets["Macrophage (positive control)"]) or "none")
print("gene overlap HSC vs GMP          :",set(sets["HSC"])&set(sets["GMP"]) or "none")

s=pd.read_csv("T4_HSPC_ssGSEA_scores.tsv",sep="\t",index_col=0)
ORDER=["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]
mono="CD14 monocyte (positive control)"; mac="Macrophage (positive control)"
X=np.column_stack([np.ones(len(s)),s[mono].values,s[mac].values])
print("\n=== after residualising on mature monocyte + macrophage content ===")
rows=[]
for m in ["HSC","LMPP","GMP","CLP","EMP (lineage control)"]:
    y=s[m].values
    beta,*_=np.linalg.lstsq(X,y,rcond=None)
    r=y-X@beta
    grp=[r[s.ecotype.values==e] for e in ORDER]
    H,p=stats.kruskal(*grp); eps=(H-2)/(len(s)-3)
    raw=stats.kruskal(*[s.loc[s.ecotype==e,m].values for e in ORDER])
    rows.append(dict(module=m,eps_raw=round((raw[0]-2)/(len(s)-3),3),KW_H_resid=round(H,2),
                     p_resid=p,eps_resid=round(eps,3),
                     **{f"resid_mean_{e}":round(g.mean(),4) for e,g in zip(ORDER,grp)}))
out=pd.DataFrame(rows); out["q_BH_resid"]=stats.false_discovery_control(out.p_resid)
print(out.to_string(index=False,float_format=lambda x:f"{x:.4g}"))
out.to_csv("T4_residualised_specificity.tsv",sep="\t",index=False)
