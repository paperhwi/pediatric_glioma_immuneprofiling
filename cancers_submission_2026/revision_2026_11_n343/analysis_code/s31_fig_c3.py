"""Redraw the ScPCA axis-linkage figure (C3) with the n = 343 reference and the new state names (figL.py logic unchanged)."""
import sys, os; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from sklearn.decomposition import PCA
A="/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_C_singlecell/ScPCA_validation/Axis_linkage/"
W=ROOT+"work_c3/"; os.makedirs(W+"deliv",exist_ok=True)
a=load(); ESo=pd.read_csv(A+"openpedcan_bulk_harmonised.tsv",sep="\t",index_col=0).loc[a.index]
imm=[c for c in ESo.columns if c not in ("ESTIMATE_immune","immune_axis","ecotype")]
mu,sd=ESo[imm].mean(),ESo[imm].std(); Zo=(ESo[imm]-mu)/sd; p=PCA(1).fit(Zo.values); ao=p.transform(Zo.values)[:,0]
if np.corrcoef(ao,ESo["T_Cell_Cytotoxicity"])[0,1]<0: ao=-ao
O=pd.DataFrame({"immune_axis":(ao-ao.mean())/ao.std(),"ecotype":a.state},index=a.index); O.to_csv(W+"openpedcan_bulk_harmonised.tsv",sep="\t")
R=pd.read_csv(ROOT+"out/T_scpca_axis_linkage_donor_tests.tsv",sep="\t"); R.to_csv(W+"link_donor_axis_state.tsv",sep="\t",index=False)
D=pd.read_csv(ROOT+"out/T_scpca_donors_projection.tsv",sep="\t",index_col=0); D.to_csv(W+"link_donors.tsv",sep="\t")
L=pd.read_csv(A+"link_donor_values_long.tsv",sep="\t"); L["axis"]=L.donor.map(D.immune_axis); L["state"]=L.donor.map(D.projected_state); L.to_csv(W+"link_donor_values_long.tsv",sep="\t",index=False)
s=open(A+"scripts/figL.py").read()
for x,y in [('COL={"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}; ST=["Immune-desert","Myeloid-dominant","Lymphocyte-inflamed"]','COL={"Immune-high":"#3B82F6","Immune-intermediate":"#EF4444","Immune-low":"#9CA3AF"}; ST=["Immune-low","Immune-intermediate","Immune-high"]'),
            ('ax.set_xticklabels(["ID","MD","LI"])','ax.set_xticklabels(["Low","Intermediate","High"])'),('"Harmonised bulk scale\\nviolin: OpenPedCan (n = 344)','"Harmonised bulk scale\\nviolin: OpenPedCan (n = 343)'),
            ('ax.text(0.02,0.98,"ID = immune-desert\\nMD = myeloid-dominant\\nLI = lymphocyte-inflamed",transform=ax.transAxes,fontsize=5,va="top"); ',''),
            ('ID<MD<LI trend P','Low<Intermediate<High trend P'),('ID n: CD8 1, CD4 0, Mye 2','Low: CD8 1, CD4 0, myeloid 2 donors'),('for s in ST: ax.scatter([],[],color=COL[s],s=14,edgecolor="k",lw=0.4,label=s)','for s in ST: ax.scatter([],[],color=COL[s],s=14,edgecolor="k",lw=0.4,label=s.replace("Immune-","").capitalize())\nax.legend(frameon=False,fontsize=5,title="Projected state",title_fontsize=5,loc="upper left")')]:
    assert x in s, x[:50]; s=s.replace(x,y)
os.chdir(W); exec(compile(s,"figL","exec"))
