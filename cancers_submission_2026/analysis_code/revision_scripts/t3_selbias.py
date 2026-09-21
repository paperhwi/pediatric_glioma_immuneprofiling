"""T3: OS-evaluable vs non-evaluable — supplementary table for Reviewer 1 item 6."""
import pandas as pd, numpy as np
from scipy import stats
ANN="/mnt/user-data/uploads/Open PBTA/Revision/FINAL MANUSCRIPT 260722 - 수정본/7. Reproducibility data/Sample annotation/sample_master_annotation.tsv"
d=pd.read_csv(ANN,sep="\t")
d["OS_eval"]=d["OS_days"].notna()&d["OS_status"].notna()
n1,n0=int(d.OS_eval.sum()),int((~d.OS_eval).sum())
rows=[]
def add(var,lev,a,b,test,stat,p):
    rows.append({"Variable":var,"Level":lev,
                 f"OS-evaluable (n={n1})":a,f"Not evaluable (n={n0})":b,
                 "Test":test,"Statistic":stat,"P":p})
for var,label in [("ecotype","Immune ecotype"),("cohort_group","Integrated molecular group"),
                  ("location_class","Anatomical location"),("reported_gender","Sex")]:
    t=pd.crosstab(d[var],d.OS_eval)
    chi2,p,dof,_=stats.chi2_contingency(t)
    first=True
    for lev in t.index:
        a,b=int(t.loc[lev,True]),int(t.loc[lev,False])
        add(label if first else "",str(lev),f"{a} ({a/n1*100:.1f}%)",f"{b} ({b/n0*100:.1f}%)",
            "Pearson chi-square" if first else "",f"chi2={chi2:.2f}, df={dof}" if first else "",
            f"{p:.4f}" if first else "")
        first=False
a=d.loc[d.OS_eval,"age_years"].dropna(); b=d.loc[~d.OS_eval,"age_years"].dropna()
u,p=stats.mannwhitneyu(a,b)
add("Age at diagnosis (years)","median [IQR]",
    f"{a.median():.1f} [{a.quantile(.25):.1f}-{a.quantile(.75):.1f}]",
    f"{b.median():.1f} [{b.quantile(.25):.1f}-{b.quantile(.75):.1f}]",
    "Mann-Whitney U",f"U={u:.0f}",f"{p:.4f}")
tab=pd.DataFrame(rows)
tab.to_csv("T3_SupplTable_OS_evaluability.tsv",sep="\t",index=False)
print(tab.to_string(index=False))

sub=[]
for g,s in d.groupby("cohort_group"):
    t=pd.crosstab(s.ecotype,s.OS_eval); chi2,p,dof,_=stats.chi2_contingency(t)
    sub.append(dict(molecular_group=g,n=len(s),chi2=round(chi2,2),df=dof,p=round(p,4),
                    **{f"%OS_eval_{k}":round(v,1) for k,v in (t[True]/t.sum(axis=1)*100).items()}))
subd=pd.DataFrame(sub); subd.to_csv("T3_within_group_ecotype_x_OSavail.tsv",sep="\t",index=False)
print("\n=== ecotype x OS-availability within molecular group ===")
print(subd.to_string(index=False))
