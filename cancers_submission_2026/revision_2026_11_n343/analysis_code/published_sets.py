"""Published gene sets replacing the ad hoc programmes (gene lists verified against the cited sources or their direct derivatives)."""
SETS={
 "T cell-inflamed GEP (Ayers 2017)":dict(genes=["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"],
     ref="Ayers M et al. J Clin Invest 2017;127:2930-40. doi:10.1172/JCI91190 (18 genes; unweighted mean of gene z-scores used here because the clinical weights are proprietary)"),
 "Cytolytic activity, CYT (Rooney 2015)":dict(genes=["GZMA","PRF1"],ref="Rooney MS et al. Cell 2015;160:48-61. doi:10.1016/j.cell.2014.12.033 (geometric mean of GZMA and PRF1)"),
 "MHC-I antigen processing and presentation, APS (Wang 2019)":dict(genes=["PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMB10","TAP1","TAP2","ERAP1","ERAP2","CANX","CALR","PDIA3","TAPBP","B2M","HLA-A","HLA-B","HLA-C"],
     ref="Wang S et al. eLife 2019;8:e49020. doi:10.7554/eLife.49020 (18 APM genes after Leone et al. 2013)"),
 "Immunogenic cell death, ICD (Garg 2016)":dict(genes=["ATG5","BAX","CALR","CASP1","CASP8","CD4","CD8A","CD8B","CXCR3","EIF2AK3","ENTPD1","FOXP3","HMGB1","HSP90AA1","IFNA1","IFNB1","IFNG","IFNGR1","IL10","IL17A","IL17RA","IL1B","IL1R1","IL6","LY96","MYD88","NLRP3","NT5E","P2RX7","PDIA3","PIK3CA","PRF1","TLR4","TNF"],
     ref="Garg AD, De Ruysscher D, Agostinis P. Oncoimmunology 2016;5:e1069938. doi:10.1080/2162402X.2015.1069938 (34 ICD-parameter genes as tabulated in e.g. Cai J et al. Front Immunol 2022;13:1011757, doi:10.3389/fimmu.2022.1011757)"),
}
# ICD core restricted to the cancer-cell-intrinsic DAMP/stress arm of the same list (subset; reported separately because the full list is dominated by immune-cell genes)
SETS["ICD, cancer-cell DAMP/ER-stress arm (subset of Garg 2016)"]=dict(genes=["CALR","HMGB1","HSP90AA1","PDIA3","EIF2AK3","ATG5","BAX","CASP8","P2RX7","ENTPD1","NT5E"],ref="subset of the Garg 2016 list (DAMP emission / ER stress / ATP metabolism members)")
CHECKPOINT_GENES=["CD274","PDCD1LG2","CD276","VTCN1","CD47","HLA-E","LGALS9","PVR","NECTIN2","TGFB1","TGFB2","IDO1","NT5E","IL10","CSF1","CCL2","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT"]
