"""Check the public n = 343 revision bundle: cohort size, state counts, exclusions, shared genes, DEG counts."""
import pandas as pd, os
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),"reproducibility_data")
a=pd.read_csv(os.path.join(D,"immune_state_assignments_n343.tsv"),sep="\t")
assert len(a)==343 and a.Kids_First_Biospecimen_ID.is_unique
assert a.state.value_counts().to_dict()=={"Immune-intermediate":144,"Immune-high":122,"Immune-low":77}
x=pd.read_csv(os.path.join(D,"T_excluded_specimens.tsv"),sep="\t"); assert len(x)==6 and not a.Kids_First_Biospecimen_ID.isin(x.Kids_First_Biospecimen_ID).any()
assert (a.state!=a.state_n344_carried).sum()==1
d=pd.read_csv(os.path.join(D,"T_adjusted_DEG_n343_all_contrasts.tsv"),sep="\t")
up=lambda c:set(d[(d.contrast==c)&(d.q_BH<0.05)&(d.adjusted_log2FC>1)].gene)
assert len(up("Immune-high vs Immune-low"))==3718 and len(up("Immune-intermediate vs Immune-low"))==431
assert len(up("Immune-high vs Immune-low")&up("Immune-intermediate vs Immune-low"))==404
assert "Kids_First_Participant_ID" not in a.columns
print("OK: n = 343; Immune-high 122 / Immune-intermediate 144 / Immune-low 77; 6 exclusions; 404 shared genes")
