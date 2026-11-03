"""Paired check of retained versus within-n343 standardized 46-feature inputs.
Both use the locked assignment protocol; k diagnostics use the separate original diagnostic protocol.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,time,itertools,shutil
import numpy as np,pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score,silhouette_score
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits

WORK=Path(__file__).parent
from analysis_paths import OUT
ORD=['Immune-high','Immune-intermediate','Immune-low']

def consensus(X,seed,B=500,k=3,n_init=10):
    rng=np.random.default_rng(seed);n=len(X);M=np.zeros((n,n));C=np.zeros((n,n))
    for b in range(B):
        idx=rng.choice(n,int(n*.8),replace=False)
        labels=KMeans(n_clusters=k,n_init=n_init,random_state=seed+b).fit_predict(X[idx])
        for ci in range(k):
            members=idx[labels==ci]
            if len(members)>=2:M[np.ix_(members,members)]+=1
        C[np.ix_(idx,idx)]+=1
    cm=np.where(C>0,M/np.maximum(C,1),0.);np.fill_diagonal(cm,1.)
    dist=1-cm;np.fill_diagonal(dist,0)
    labels=fcluster(linkage(squareform(dist,checks=False),'average'),t=k,criterion='maxclust')
    return labels,cm,dist

def job(kind,seed_or_k):
    threadpool_limits(limits=1)
    F=pd.read_csv(OUT/'inputs/retained_z_n343.tsv',sep='\t',index_col=0)
    Z=pd.read_csv(OUT/'inputs/restandardized_z_n343.tsv',sep='\t',index_col=0)
    start=time.time();(OUT/'.cache/seeds').mkdir(parents=True,exist_ok=True)
    if kind in ['retained','restandardized']:
        X=F if kind=='retained' else Z
        labels,cm,dist=consensus(X.values,seed_or_k)
        means=X.assign(cluster=labels).groupby('cluster')['ssGSEA_T_Cell_Cytotoxicity'].mean().sort_values()
        names=dict(zip(means.index,['Immune-low','Immune-intermediate','Immune-high']))
        named=pd.Series(labels,index=X.index).map(names)
        named.to_frame(str(seed_or_k)).to_csv(OUT/'.cache/seeds'/f'seed_{kind}_{seed_or_k}.tsv',sep='\t')
        if seed_or_k==42:np.savez_compressed(OUT/'tables'/f'consensus_{kind}_seed42.npz',consensus=cm,ids=np.array(X.index,dtype=str))
        return {'kind':kind,'seed':seed_or_k,'seconds':time.time()-start,'counts':named.value_counts().to_dict()}
    if kind=='ksweep':
        k=seed_or_k;labels,cm,dist=consensus(Z.values,42,B=1000,k=k,n_init=1)
        up=np.sort(cm[np.triu_indices(len(Z),1)]);cdf=np.arange(1,len(up)+1)/len(up)
        sizes=pd.Series(labels).value_counts();same=np.equal.outer(labels,labels)&~np.eye(len(Z),dtype=bool)
        return {'kind':kind,'k':k,'PAC':float(((up>.1)&(up<.9)).mean()),'silhouette_consensus':float(silhouette_score(dist,labels,metric='precomputed')),
                'consensus_CDF_AUC':float(np.sum(np.diff(np.r_[0.,up])*cdf)),'within_consensus':float(cm[same].mean()),
                'min_cluster_size':int(sizes.min()),'sizes':';'.join(map(str,sorted(sizes,reverse=True))),'seconds':time.time()-start}
    if kind=='gap':
        X=Z.values;n,p=X.shape;ks=list(range(1,11));rng=np.random.default_rng(42)
        pca=PCA().fit(X);xp=X@pca.components_.T;lo,hi=xp.min(0),xp.max(0)
        def wk(D,k):
            return float(((D-D.mean(0))**2).sum()) if k==1 else float(KMeans(n_clusters=k,n_init=10,random_state=42).fit(D).inertia_)
        logw=np.log([wk(X,k) for k in ks]);ref=np.zeros((50,len(ks)))
        for b in range(50):
            d=rng.uniform(lo,hi,size=(n,p))@pca.components_;ref[b]=np.log([wk(d,k) for k in ks])
        g=ref.mean(0)-logw;sk=ref.std(0)*np.sqrt(1+1/50)
        table=pd.DataFrame({'k':ks,'logW':logw,'gap':g,'s_k':sk,'meets_1SE_rule':[bool(g[i]>=g[i+1]-sk[i+1]) for i in range(9)]+[False]})
        table.to_csv(OUT/'tables/gap_restandardized_n343.tsv',sep='\t',index=False)
        selected=table.loc[table.meets_1SE_rule,'k']
        return {'kind':'gap','selected_k':int(selected.iloc[0]) if len(selected) else None,'seconds':time.time()-start}

def main():
    start=time.time();F=pd.read_csv(OUT/'inputs/retained_z_n343.tsv',sep='\t',index_col=0);Z=pd.read_csv(OUT/'inputs/restandardized_z_n343.tsv',sep='\t',index_col=0)
    obs=pd.read_csv(OUT/'inputs/public_annotations_n343.tsv',sep='\t',index_col=0)
    plan={'cluster_assignment':{'features':46,'n':343,'seeds':list(range(42,67)),'B':500,'pItem':.8,'KMeans_n_init':10,'linkage':'average','distance':'1-consensus','k':3,'naming':'ordered mean T-cell cytotoxicity','aggregate':'per-sample mode; pandas lexical tie rule'},
          'rescaling':{'formula':'(retained_Z - mean_n343) / sample_SD_n343','ddof':1,'raw_affine_equivalence':'re-z-scoring an affine standardized feature equals standardizing its original score, up to floating precision; no re-scoring of genes'},
          'k_diagnostics':{'ks':list(range(2,11)),'B':1000,'pItem':.8,'KMeans_n_init':1,'seed':42,'gap_references':50,'gap_KMeans_n_init':10,'gap_rule':'first k satisfying 1SE'},
          'comparison':'paired same seeds and sample order; reproduce old modal assignments before interpreting scaling sensitivity'}
    (OUT/'code/analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    scale=pd.DataFrame({'feature':F.columns,'retained_mean_n343':F.mean().values,'retained_SD_n343':F.std().values,'new_mean_n343':Z.mean().values,'new_SD_n343':Z.std().values})
    scale.to_csv(OUT/'tables/feature_standardization_audit.tsv',sep='\t',index=False)
    tasks=[(kind,seed) for kind in ['retained','restandardized'] for seed in range(42,67)]+[('ksweep',k) for k in range(2,11)]+[('gap',0)]
    records=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(job,*task):task for task in tasks}
        for future in as_completed(futures):
            result=future.result();records.append(result);print(json.dumps(result),flush=True)
            (WORK/'cluster_progress.json').write_text(json.dumps({'completed':len(records),'total':len(tasks),'elapsed_s':time.time()-start,'latest':result},indent=2),encoding='utf-8')
    seedframes={};modes={};summary=[]
    for kind in ['retained','restandardized']:
        S=pd.concat([pd.read_csv(OUT/'.cache/seeds'/f'seed_{kind}_{seed}.tsv',sep='\t',index_col=0) for seed in range(42,67)],axis=1).loc[obs.index]
        S.to_csv(OUT/'tables'/f'assignments_25seeds_{kind}.tsv',sep='\t')
        mode=S.mode(axis=1)[0];freq=S.eq(mode,axis=0).mean(axis=1);ari=[adjusted_rand_score(S[str(i)],S[str(j)]) for i,j in itertools.combinations(range(42,67),2)]
        ties=S.apply(lambda row:len(row.mode()),axis=1)>1
        summary.append({'input':kind,'n':343,**{f'n_{s}':int((mode==s).sum()) for s in ORD},'median_pairwise_ARI':float(np.median(ari)),'min_pairwise_ARI':float(min(ari)),'max_pairwise_ARI':float(max(ari)),
                        'same_label_all25':int((freq==1).sum()),'min_modal_frequency':float(freq.min()),'modal_ties':int(ties.sum()),'agreement_vs_delivered':float((mode==obs.immune_state).mean()),'ARI_vs_delivered':float(adjusted_rand_score(mode,obs.immune_state))})
        seedframes[kind]=S;modes[kind]=mode
        obs[f'{kind}_modal']=mode;obs[f'{kind}_modal_frequency']=freq
    reproduced=bool((modes['retained']==obs.immune_state).all())
    if not reproduced:raise RuntimeError('Delivered old modal labels did not reproduce: report platform/procedure difference before interpreting scaling comparison.')
    summary=pd.DataFrame(summary);summary.to_csv(OUT/'tables/clustering_summary.tsv',sep='\t',index=False)
    ct=pd.crosstab(obs.immune_state.rename('delivered_state'),modes['restandardized'].rename('restandardized_state')).reindex(index=ORD,columns=ORD,fill_value=0)
    ct.to_csv(OUT/'tables/state_transition_counts.tsv',sep='\t')
    obs['assignment_changed']=modes['restandardized']!=obs.immune_state
    obs[['immune_state','retained_modal','retained_modal_frequency','restandardized_modal','restandardized_modal_frequency','assignment_changed','molecular_group','RNA_library']].to_csv(OUT/'tables/per_sample_clustering_comparison.tsv',sep='\t')
    obs.loc[obs.assignment_changed,['immune_state','restandardized_modal','restandardized_modal_frequency','molecular_group','RNA_library']].to_csv(OUT/'tables/changed_assignments.tsv',sep='\t')
    rows=[]
    for seed in range(42,67):
        a=seedframes['retained'][str(seed)];b=seedframes['restandardized'][str(seed)]
        rows.append({'seed':seed,'paired_ARI':adjusted_rand_score(a,b),'paired_agreement':float((a==b).mean()),'retained_ARI_vs_delivered':adjusted_rand_score(a,obs.immune_state),'new_ARI_vs_new_modal':adjusted_rand_score(b,modes['restandardized']),**{f'new_n_{s}':int((b==s).sum()) for s in ORD}})
    pd.DataFrame(rows).to_csv(OUT/'tables/seed_comparison.tsv',sep='\t',index=False)
    sweep=pd.DataFrame([r for r in records if r['kind']=='ksweep']).sort_values('k').drop(columns=['kind','seconds']);sweep['delta_AUC']=sweep.consensus_CDF_AUC.pct_change()
    sweep.to_csv(OUT/'tables/ksweep_restandardized_n343.tsv',sep='\t',index=False)
    result={'baseline_modal_reproduced':reproduced,'modal_agreement':float((modes['restandardized']==obs.immune_state).mean()),'modal_ARI':adjusted_rand_score(modes['restandardized'],obs.immune_state),'changed_samples':int(obs.assignment_changed.sum()),
            'best_PAC_k':int(sweep.loc[sweep.PAC.idxmin(),'k']),'best_silhouette_k':int(sweep.loc[sweep.silhouette_consensus.idxmax(),'k']),'gap_selected_k':next(r['selected_k'] for r in records if r['kind']=='gap'),'elapsed_s':time.time()-start,'clustering_summary':summary.to_dict('records')}
    (OUT/'tables/clustering_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('FINAL',json.dumps(result),flush=True)

if __name__=='__main__':main()
