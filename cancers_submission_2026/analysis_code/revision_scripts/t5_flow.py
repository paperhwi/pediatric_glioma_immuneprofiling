"""Figure 1B — analysis-set derivation. Branching, because the location and survival
subsets are not nested (only 180 of the 251 survival-evaluable tumours are in the n=258 set)."""
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42})

NAVY,BLUE,RED,GREY="#1E3A5F","#1D4ED8","#B91C1C","#6B7280"
fig,ax=plt.subplots(figsize=(9.2,7.0),dpi=300)
ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")

def node(cx,cy,w,h,lines,fc="#FFFFFF",ec=NAVY,lw=1.15):
    """lines = [(text, fontsize, weight, color), ...] stacked and centred."""
    ax.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,boxstyle="round,pad=0.55,rounding_size=1.1",
                                fc=fc,ec=ec,lw=lw,zorder=2))
    tot=sum(l[1] for l in lines); gap=h/(tot+len(lines)*0.9)
    y=cy+h/2-gap*1.5
    for txt,fs,wt,col in lines:
        ax.text(cx,y,txt,ha="center",va="center",fontsize=fs,weight=wt,color=col,zorder=3,linespacing=1.35)
        y-=gap*(fs*0.92)
    return cy-h/2, cy+h/2

def arrow(x1,y1,x2,y2,color=NAVY,lw=1.1):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,
                                 color=color,lw=lw,shrinkA=0,shrinkB=0,zorder=1))

def excl(cx,cy,w,h,txt):
    ax.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,boxstyle="round,pad=0.45,rounding_size=0.9",
                                fc="#FBF3F3",ec=RED,lw=0.95,linestyle=(0,(3,2)),zorder=2))
    ax.text(cx,cy,txt,ha="center",va="center",fontsize=7.0,color="#7F1D1D",zorder=3,linespacing=1.4)

# ---- source
node(50,93,44,8,[("OpenPedCan v15",9.0,"bold",NAVY),
                 ("paediatric CNS tumour biospecimens",7.6,"normal","#374151")],fc="#EDF2F9")
arrow(50,89,50,84.5)
excl(84,88.6,28,8.2,"Excluded: non-primary and repeat\nbiospecimens; diagnoses outside the\nfour integrated groups")
ax.add_patch(FancyArrowPatch((50.8,88.6),(70,88.6),arrowstyle="-|>",mutation_scale=11,color=RED,lw=0.95,
                             shrinkA=0,shrinkB=0,zorder=1))

# ---- main cohort
node(50,78,58,10,[("Main cohort   n = 349",10.0,"bold",NAVY),
                  ("one independent primary biospecimen per patient",7.6,"normal","#374151"),
                  ("DMG_K27 175  ·  pHGG_WT 125  ·  DHG_G34 31  ·  IHG 18",7.4,"normal","#374151")],
     fc="#E4EBF6")
ax.text(50,70.6,"immune deconvolution, ssGSEA scoring and consensus clustering used all 349 tumours",
        ha="center",va="center",fontsize=7.4,style="italic",color="#4B5563",zorder=4,
        bbox=dict(boxstyle="round,pad=0.25",fc="white",ec="none"))

# ---- split
ax.plot([50,50],[73,68.2],color=NAVY,lw=1.1,zorder=1)
ax.plot([27,73],[64.4,64.4],color=NAVY,lw=1.1,zorder=1)
ax.plot([50,50],[68.2,64.4],color=NAVY,lw=1.1,zorder=1)
arrow(27,64.4,27,59.6); arrow(73,64.4,73,59.6)

# ---- left branch
node(27,54,40,10,[("Anatomical-location analyses",8.4,"bold",NAVY),
                  ("n = 332",9.8,"bold",BLUE),
                  ("contingency tests and Cramér's V",7.3,"normal","#374151")])
excl(5,54,12,8,"− 17\nno usable\nlocation")
ax.add_patch(FancyArrowPatch((7,54),(11,54),arrowstyle="-|>",mutation_scale=11,color=RED,lw=0.95,shrinkA=0,shrinkB=0,zorder=1))
arrow(27,49,27,43.4)
node(27,38,40,10,[("Sequential PERMANOVA",8.4,"bold",NAVY),
                  ("n = 258",9.8,"bold",BLUE),
                  ("midline, hemispheric, posterior fossa",7.3,"normal","#374151")])
excl(5,38,12,8,"− 74\nmulti-\ncompartment")
ax.add_patch(FancyArrowPatch((7,38),(11,38),arrowstyle="-|>",mutation_scale=11,color=RED,lw=0.95,shrinkA=0,shrinkB=0,zorder=1))

# ---- right branch
node(73,54,40,10,[("Survival analyses",8.4,"bold",NAVY),
                  ("n = 251   (208 deaths)",9.8,"bold",BLUE),
                  ("Kaplan–Meier and multivariable Cox",7.3,"normal","#374151")])
excl(95,54,10,8,"− 98\nno OS\nannotation")
ax.add_patch(FancyArrowPatch((93,54),(90,54),arrowstyle="<|-",mutation_scale=11,color=RED,lw=0.95,shrinkA=0,shrinkB=0,zorder=1))
arrow(73,49,73,43.4)
node(73,38,40,10,[("Missingness sensitivity",8.4,"bold",NAVY),
                  ("n = 349",9.8,"bold",BLUE),
                  ("inverse-probability-of-observation weighting",7.1,"normal","#374151")])

# ---- non-nesting note
ax.plot([27,27],[33,28.5],color=GREY,lw=0.9,ls=":",zorder=1)
ax.plot([73,73],[33,28.5],color=GREY,lw=0.9,ls=":",zorder=1)
ax.plot([27,73],[28.5,28.5],color=GREY,lw=0.9,ls=":",zorder=1)
ax.plot([50,50],[28.5,25.5],color=GREY,lw=0.9,ls=":",zorder=1)
ax.add_patch(FancyBboxPatch((7,7),86,18,boxstyle="round,pad=0.6,rounding_size=1.2",
                            fc="#FFF8E7",ec="#C08A2E",lw=1.0,zorder=2))
ax.text(50,21,"The two branches are not nested",ha="center",va="center",
        fontsize=8.4,weight="bold",color="#7A5512",zorder=3)
ax.text(50,13.6,"Of the 251 survival-evaluable tumours, 180 fall within the n = 258 location subset and 3 have no usable location\n"
                "annotation. The location and survival analyses therefore draw different subsets of the same 349 tumours and\n"
                "must not be read as a single attrition sequence. Supplementary Table S17 reconciles the counts.",
        ha="center",va="center",fontsize=7.3,color="#5B4210",zorder=3,linespacing=1.6)

for e in ("png","pdf"):
    fig.savefig(f"Figure1B_attrition.{e}",dpi=300,bbox_inches="tight",facecolor="white")
print("saved")
