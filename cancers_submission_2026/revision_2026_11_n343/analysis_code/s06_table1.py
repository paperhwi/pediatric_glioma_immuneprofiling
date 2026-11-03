"""Gap G2 — revised Table 1 (Reviewer 1: duplicated count rows; title says "molecular" but no molecular rows).
n = 343 (single exclusion rule). Columns: all tumours and each immune state.
Clinical rows and molecular rows (DNA-based calls from OpenPedCan v15 MAF/CNV, as in Table S28); every molecular row reports
k/n with its own denominator (tumours with the relevant DNA data). P: Kruskal-Wallis (age) or Monte-Carlo chi-square
(fixed margins, 9,999 permutations of the state labels, seed 20260722); BH q across molecular rows."""
import sys, json, os; sys.path.insert(0, "/home/claude/q343/src")
from common import *
from scipy import stats
from lifelines import KaplanMeierFitter
O = ROOT+"out/"
P = load(); M = pd.read_pickle(ROOT+"markers343.pkl").loc[P.index]
st = P.state; COLS = ["All"] + ORD; rng = np.random.default_rng(20260722)
def _chi2(ci, gi, nc, ng):
    tab = np.bincount(ci * ng + gi, minlength=nc * ng).reshape(nc, ng).astype(float)
    e = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum(); m = e > 0
    return (((tab - e) ** 2)[m] / e[m]).sum()
def mc_chi2(cat, grp, B=9999):
    ok = cat.notna() & grp.notna(); c = pd.factorize(cat[ok].astype(str))[0]; g = pd.factorize(grp[ok])[0]
    if len(set(c)) < 2: return np.nan
    nc, ng = c.max() + 1, g.max() + 1; obs = _chi2(c, g, nc, ng)
    cnt = sum(_chi2(c, rng.permutation(g), nc, ng) >= obs - 1e-9 for _ in range(B))
    return (cnt + 1) / (B + 1)
ROWS = []
def add(section, label, value=None, p=np.nan, test="", note=""):
    ROWS.append(dict(section=section, characteristic=label, **(value or {c: "" for c in COLS}), P=p, test=test, note=note))
def cat_rows(section, label, s, order=None, note=""):
    s = s.copy(); p = mc_chi2(s, st); order = order or list(s.dropna().value_counts().index)
    add(section, label, None, p, "Monte-Carlo χ²", note)
    for lv in order:
        v = {}
        for c in COLS:
            ss = s if c == "All" else s[st == c]; n = ss.notna().sum(); k = (ss == lv).sum()
            v[c] = f"{k} ({100 * k / n:.1f}%)" if n else "–"
        add(section, f"  {lv}", v)
    nm = s.isna().sum()
    if nm: add(section, "  Not available", {c: str(int((s if c == "All" else s[st == c]).isna().sum())) for c in COLS})
def bin_row(section, label, s, note=""):
    p = mc_chi2(s, st); v = {}
    for c in COLS:
        ss = s if c == "All" else s[st == c]; n = int(ss.notna().sum()); k = int((ss == 1).sum())
        v[c] = f"{k}/{n} ({100 * k / n:.1f}%)" if n else "–"
    add(section, label, v, p, "Monte-Carlo χ²", note)
# ---- header counts (one row only: n per column, replaces the two duplicated rows)
add("", "Tumours (one independent primary RNA specimen per patient), n", {c: str(int(len(P) if c == "All" else (st == c).sum())) for c in COLS})
# ---- clinical
age = P.age_first_dx_years; p_age = stats.kruskal(*[age[st == e].dropna() for e in ORD]).pvalue
add("Clinical", "Age at first diagnosis, years, median (IQR)", {c: (lambda a: f"{a.median():.1f} ({a.quantile(.25):.1f}–{a.quantile(.75):.1f})")(age if c == "All" else age[st == c]) for c in COLS}, p_age, "Kruskal–Wallis")
ag = pd.cut(age, [-0.01, 3, 10, 15, 40], labels=["< 3", "3 to < 10", "10 to < 15", "≥ 15"]).astype(object)
cat_rows("Clinical", "Age group, years", ag, ["< 3", "3 to < 10", "10 to < 15", "≥ 15"])
sex = P.reported_gender.where(P.reported_gender.isin(["Female", "Male"])); cat_rows("Clinical", "Sex", sex, ["Female", "Male"])
td = P.tumor_descriptor.replace({"Primary Tumor": "Initial CNS Tumor", "Deceased": "Autopsy (post-mortem)", "Second Malignancy": "Recurrence/progressive",
                                  "Recurrence": "Recurrence/progressive", "Progressive": "Recurrence/progressive"})
cat_rows("Clinical", "Specimen timing", td, ["Initial CNS Tumor", "Recurrence/progressive", "Autopsy (post-mortem)"])
reg = P.CNS_region.replace({"Ventricles": "Other", "Suprasellar": "Other", "Spine": "Spine"})
cat_rows("Clinical", "Anatomical region (OpenPedCan CNS_region)", reg, ["Midline", "Hemispheric", "Mixed", "Posterior fossa", "Spine", "Other"])
eor = P.extent_of_tumor_resection.replace({"Unavailable": np.nan}); cat_rows("Clinical", "Extent of resection", eor, ["Gross/Near total resection", "Partial resection", "Biopsy only"])
lib = P.RNA_library; cat_rows("Clinical", "RNA library", lib, ["stranded", "poly-A stranded", "poly-A", "exome_capture"])
# survival
v = {}
for c in COLS:
    d = P[P.os_eval] if c == "All" else P[P.os_eval & (st == c)]
    km = KaplanMeierFitter().fit(d.time / 365.25, d.event); v[c] = f"{len(d)} / {int(d.event.sum())} / {km.median_survival_time_:.2f}"
add("Clinical", "Overall survival: evaluable n / deaths / median OS, years (corrected OS)", v, np.nan, "", "Log-rank and Cox results: Figure 7")
# ---- molecular
add("Molecular", "Integrated molecular group (WHO CNS5)", None, mc_chi2(P.cohort_group, st), "Monte-Carlo χ² (n = 343)", "")
for lv, labl in [("DMG_K27", "Diffuse midline glioma, H3 K27-altered"), ("DHG_G34", "Diffuse hemispheric glioma, H3 G34-mutant"),
                 ("pHGG_WT", "Diffuse pediatric-type HGG, H3-wildtype and IDH-wildtype"), ("IHG", "Infant-type hemispheric glioma")]:
    add("Molecular", f"  {labl}", {c: (lambda s: f"{(s == lv).sum()} ({100 * (s == lv).mean():.1f}%)")(P.cohort_group if c == "All" else P.cohort_group[st == c]) for c in COLS})
dna = M.has_DNA.astype(int); add("Molecular", "Tumours with matched DNA (WGS/WXS), n", {c: str(int((dna if c == "All" else dna[st == c]).sum())) for c in COLS})
MOL = [("TP53_alt", "TP53 altered (mutation or deep deletion)"), ("ATRX_mut", "ATRX mutation"), ("PPM1D_mut", "PPM1D mutation"), ("ACVR1_mut", "ACVR1 mutation"),
       ("PI3K_pathway_alt", "PI3K pathway (PIK3CA, PIK3R1, PTEN)"), ("NF1_alt", "NF1 altered"), ("BCOR_mut", "BCOR mutation"), ("FGFR1_mut", "FGFR1 mutation"),
       ("BRAF_mut", "BRAF mutation"), ("SETD2_mut", "SETD2 mutation"), ("PDGFRA_alt", "PDGFRA altered (mutation or amplification)"), ("EGFR_amp", "EGFR amplification"),
       ("MYCN_amp", "MYCN amplification"), ("MYC_amp", "MYC amplification"), ("CDK4_amp", "CDK4 amplification"), ("CDK6_amp", "CDK6 amplification"), ("CCND2_amp", "CCND2 amplification"),
       ("MDM4_amp", "MDM4 amplification"), ("CDKN2A/B_deep_del", "CDKN2A/B homozygous deletion"), ("RB1_alt", "RB1 altered"), ("MMR/POLE_mut", "MMR or POLE mutation"),
       ("Hypermutant_TMB>=10", "Hypermutant (TMB ≥ 10 mut/Mb)"), ("MGMT_methylated", "MGMT promoter methylated (methylation array)"), ("Cancer_predisposition", "Reported cancer predisposition")]
mol_idx = []
for col, labl in MOL:
    s = M[col].copy()
    if col == "Cancer_predisposition": s = s.astype(float)
    if s.dropna().sum() == 0: continue
    mol_idx.append(len(ROWS)); bin_row("Molecular", labl, s)
T = pd.DataFrame(ROWS); T["q_BH_molecular"] = np.nan
pm = T.loc[mol_idx, "P"]; T.loc[mol_idx, "q_BH_molecular"] = bh(pm.values)
T = T.rename(columns={e: f"{e} (n = {(st == e).sum()})" for e in ORD} | {"All": f"All (n = {len(P)})"})
T.to_csv(O + "T_Table1_n343.tsv", sep="\t", index=False)
S = dict(title="Table 1. Clinical and molecular characteristics of the analysis cohort (n = 343), overall and by immune state",
         significant_q05=T.loc[T.q_BH_molecular < 0.05, ["characteristic", "P", "q_BH_molecular"]].to_dict("records"),
         nominal_p05=T.loc[(T.P < 0.05), ["characteristic", "P"]].to_dict("records"))
json.dump(S, open(O + "T_Table1_summary.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250); print(T.drop(columns=["section", "note"]).to_string()); print(json.dumps(S, indent=1, default=float))
