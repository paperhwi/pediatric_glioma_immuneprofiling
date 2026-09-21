"""Audit: two signatures sit in themes that contradict their gene content and their
single-cell attribution. Quantify the impact on the published theme-composite result."""
import numpy as np, pandas as pd
from scipy import stats
ANN="/mnt/user-data/uploads/Open PBTA/Revision/FINAL MANUSCRIPT 260722 - 수정본/7. Reproducibility data/Sample annotation/sample_master_annotation.tsv"
d=pd.read_csv(ANN,sep="\t").set_index("sample")
sig=[c for c in d.columns if c.startswith("ssGSEA_")]
Z=(d[sig]-d[sig].mean())/d[sig].std()
eco=d["ecotype"]; ORDER=["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]

PUBLISHED={
 "T-cell axis":["T_Cell_Cytotoxicity","T_Cell_Exhaustion","Tregs_Friebel2020","Chemokine_T_Cell_Recruitment",
                "NK_Cell_Activity","Dendritic_Cell_Activation"],
 "Antigen presentation":["MHC_Class_I","MHC_Class_II"],
 "IFN / chemokine":["IFN_Gamma_Response","IFN_Alpha_Response"],
 "Myeloid / microglia":["Microglia_Core_Homeostatic","Microglia_Klemm2020","MDM_Klemm2020","MgTAM_Antunes2021",
                        "MoTAM_Antunes2021","DAM_KerenShaul2017","M1_Macrophage","M2_Macrophage",
                        "Neutrophil_Activation","Glioma_Inflammatory_Wang2017"],
 "Signaling / suppression":["MAPK_Activity","TGFb_Immunosuppression","Cell_Cycle_Proliferation","Stemness_Brain_Tumor"],
}
CORRECTED={k:list(v) for k,v in PUBLISHED.items()}
CORRECTED["Myeloid / microglia"].remove("Glioma_Inflammatory_Wang2017")
CORRECTED["Myeloid / microglia"].append("Dendritic_Cell_Activation")
CORRECTED["T-cell axis"].remove("Dendritic_Cell_Activation")
CORRECTED["T-cell axis"].append("Glioma_Inflammatory_Wang2017")

def composites(themes):
    return pd.DataFrame({t:Z[[f"ssGSEA_{s}" for s in v]].mean(axis=1) for t,v in themes.items()})
def kw(T):
    out=[]
    for t in T.columns:
        g=[T.loc[eco==e,t].values for e in ORDER]
        H,p=stats.kruskal(*g); eps=(H-2)/(len(T)-3)
        out.append(dict(theme=t,epsilon_sq=round(eps,3),p=p,
                        **{e.split("-")[0][:4]:round(float(np.mean(x)),3) for e,x in zip(ORDER,g)}))
    return pd.DataFrame(out)

P=kw(composites(PUBLISHED)); C=kw(composites(CORRECTED))
m=P.merge(C,on="theme",suffixes=("_published","_corrected"))
print("=== theme composite, Kruskal-Wallis across ecotypes ===")
print(m[["theme","epsilon_sq_published","epsilon_sq_corrected"]].to_string(index=False))
print("\n=== published ranking ==="); print(P.sort_values("epsilon_sq",ascending=False)[["theme","epsilon_sq"]].to_string(index=False))
print("\n=== corrected ranking ==="); print(C.sort_values("epsilon_sq",ascending=False)[["theme","epsilon_sq"]].to_string(index=False))

# does the Myeloid-dominant label survive? myeloid theme rank of each ecotype
Tp=composites(PUBLISHED); Tc=composites(CORRECTED)
print("\n=== mean myeloid/microglia composite by ecotype ===")
for nm,T in [("published",Tp),("corrected",Tc)]:
    v=T.groupby(eco)["Myeloid / microglia"].mean().reindex(ORDER).round(3)
    print(f"  {nm:10} " + "  ".join(f"{k}={x}" for k,x in v.items()))
print("\n=== mean T-cell axis composite by ecotype ===")
for nm,T in [("published",Tp),("corrected",Tc)]:
    v=T.groupby(eco)["T-cell axis"].mean().reindex(ORDER).round(3)
    print(f"  {nm:10} " + "  ".join(f"{k}={x}" for k,x in v.items()))
# correlation between the two themes, before and after
print("\n=== myeloid x T-cell theme correlation (the contamination effect) ===")
for nm,T in [("published",Tp),("corrected",Tc)]:
    r=np.corrcoef(T["Myeloid / microglia"],T["T-cell axis"])[0,1]
    print(f"  {nm:10} Pearson r = {r:.3f}")
m.to_csv("W2_theme_audit.tsv",sep="\t",index=False)
