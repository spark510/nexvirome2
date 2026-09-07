"""Nonnegative multi-view factorization with explicitly masked missing entries."""
from dataclasses import dataclass
from collections import defaultdict
import math
from .common import rows,stage,table,dump,sha


def propose(settings,model,coordinates,output):
    from pathlib import Path
    from .common import config
    model=Path(model)
    with stage(output,settings,[model/'loadings.tsv',model/'diagnostics.json',coordinates]) as (out,_):
        diagnostics=config(model/'diagnostics.json')
        source=rows(model/'loadings.tsv')
        loadings={r['entity']:r for r in source}
        if len(loadings)!=len(source):
            raise ValueError('Duplicate latent entity')
        by_contig=defaultdict(list)
        inactive=set(diagnostics.get('inactive_entities',[]))
        if inactive-set(loadings):
            raise ValueError('Unknown inactive latent entity')
        seen=set()
        for row in rows(coordinates):
            entity=row['entity']
            if entity not in loadings or entity in seen:
                raise ValueError('Unknown/duplicate coordinate entity')
            seen.add(entity)
            start,end=int(row['start']),int(row['end'])
            if not row['contig'] or not 0<=start<end:
                raise ValueError('Invalid original coordinates')
            vector=[float(v) for k,v in loadings[entity].items() if k.startswith('component_')]
            is_inactive=entity in inactive
            if not vector or any(not math.isfinite(v) or not 0<=v<=1 for v in vector):
                raise ValueError('Invalid normalized loadings')
            if is_inactive:
                if any(vector) or loadings[entity].get('status')!='inactive':
                    raise ValueError('Invalid inactive loadings')
            elif abs(sum(vector)-1)>1e-6 or loadings[entity].get('status','active')!='active':
                raise ValueError('Invalid normalized loadings')
            by_contig[row['contig']].append((start,end,None if is_inactive else vector))
        if seen!=set(loadings):
            raise ValueError('Coordinates must cover every latent entity')
        result=[]
        minimum=float(settings.get('min_loading',.7))
        if not 0<=minimum<=1:
            raise ValueError('min_loading must be a fraction')
        for contig,windows in by_contig.items():
            windows.sort()
            for left,right in zip(windows,windows[1:]):
                if left[1]>right[0]:
                    raise ValueError('Overlapping latent coordinate windows')
                if not diagnostics.get('proposal_eligible') or left[1]!=right[0]:
                    continue
                a,b=left[2],right[2]
                if a is None or b is None:
                    continue
                if a.index(max(a))!=b.index(max(b)) and min(max(a),max(b))>=minimum:
                    result.append(dict(contig=contig,start=left[1],end=right[0],method='multiview',
                         score=min(max(a),max(b))*float(diagnostics['stability']),
                         detail='Stable mixed-feature switch; requires independent read/graph support'))
        table(out/'proposals.tsv',result,['contig','start','end','method','score','detail'])


@dataclass
class FeatureMatrix:
    entities: list
    columns: list
    values: object
    observed: object

    @classmethod
    def load(cls,path):
        import numpy as np
        records=rows(path)
        entities=sorted({r['entity'] for r in records})
        columns=sorted({(r['block'],r['feature']) for r in records})
        if len(entities)<2 or not columns or any(not x for x in entities) or any(not all(c) for c in columns):
            raise ValueError('Need at least two entities and named feature blocks')
        ei,ci={v:i for i,v in enumerate(entities)},{v:i for i,v in enumerate(columns)}
        values=np.zeros((len(entities),len(columns)))
        observed=np.zeros_like(values,dtype=bool)
        seen=set()
        for r in records:
            key=(ei[r['entity']],ci[(r['block'],r['feature'])])
            if key in seen:
                raise ValueError('Duplicate entity/block/feature')
            seen.add(key)
            if r['value'] not in ('','NA'):
                value=float(r['value'])
                if not math.isfinite(value) or value<0:
                    raise ValueError('Features must be finite and nonnegative')
                values[key]=value
                observed[key]=True
        if not observed.any(axis=1).all():
            raise ValueError('Entities without any measured feature cannot be factorized')
        return cls(entities,columns,values,observed)


def normalize(matrix,settings):
    import numpy as np
    x=np.zeros_like(matrix.values)
    weights=np.zeros_like(x)
    scales=[]
    blocks=defaultdict(list)
    for j,(block,feature) in enumerate(matrix.columns):
        blocks[block].append(j)
        observed=matrix.values[matrix.observed[:,j],j]
        active=len(observed)>=2 and np.ptp(observed)>1e-12
        scale=float(np.quantile(observed,.95)) if active else 1.0
        x[:,j]=np.clip(matrix.values[:,j]/scale,0,1) if active else 0
        scales.append(dict(block=block,feature=feature,scale=scale,active=bool(active)))
    configured=settings.get('block_weights',{})
    if set(configured)-set(blocks):
        raise ValueError('Unknown configured feature block')
    for block,columns in blocks.items():
        weight=float(configured.get(block,1))
        if not math.isfinite(weight) or weight<0:
            raise ValueError('Block weights must be finite and nonnegative')
        active=[j for j in columns if scales[j]['active']]
        if active:
            present=matrix.observed[:,active]
            denominator=np.maximum(1,present.sum(axis=1))
            weights[:,active]=present*weight/denominator[:,None]
    if not (weights.sum(axis=1)>0).all():
        raise ValueError('An entity has no observed variable feature with positive block weight')
    return x,weights,scales


def factorize(x,weights,rank,seed,max_iter=2000,tolerance=1e-6):
    import numpy as np
    if rank<2 or rank>min(x.shape) or max_iter<2 or not 0<tolerance<1:
        raise ValueError('Invalid NMF rank/iteration/tolerance')
    rng=np.random.default_rng(seed)
    w=rng.uniform(.1,1,(x.shape[0],rank))
    h=rng.uniform(.1,1,(rank,x.shape[1]))
    previous=float('inf')
    for iteration in range(max_iter):
        h*= (w.T@(weights*x))/(w.T@(weights*(w@h))+1e-12)
        w*= ((weights*x)@h.T)/((weights*(w@h))@h.T+1e-12)
        scale=np.maximum(h.sum(axis=1),1e-12)
        h/=scale[:,None]
        w*=scale[None,:]
        loss=float(np.sum(weights*(x-w@h)**2))
        if not np.isfinite(loss):
            raise ValueError('NMF diverged')
        if abs(previous-loss)<=tolerance*max(1e-12,previous) and math.isfinite(previous):
            break
        previous=loss
    return w,h,dict(loss=loss,iterations=iteration+1,converged=iteration+1<max_iter)


def fit(settings,features,output):
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    with stage(output,settings,[features]) as (out,_):
        matrix=FeatureMatrix.load(features)
        x,weights,scales=normalize(matrix,settings)
        active=(weights*x).sum(axis=1)>0
        seeds=settings.get('seeds',[17,29,43])
        if len(set(seeds))<2:
            raise ValueError('Multiple initialization seeds are required for stability assessment')
        rank=int(settings.get('rank',2))
        active_columns=sum(s['active'] for s in scales)
        if rank>min(int(active.sum()),active_columns):
            raise ValueError('Rank exceeds entities or variable features')
        fits=[factorize(x[active],weights[active],rank,int(seed),int(settings.get('max_iter',2000)),
                        float(settings.get('tolerance',1e-6))) for seed in seeds]
        best=min(range(len(fits)),key=lambda i:(fits[i][2]['loss'],i))
        w,h,diagnostics=fits[best]
        similarities=[]
        for i,(other,_,_) in enumerate(fits):
            if i==best:
                continue
            scores=(w.T@other)/(np.maximum(np.linalg.norm(w,axis=0),1e-12)[:,None]*
                                np.maximum(np.linalg.norm(other,axis=0),1e-12)[None,:])
            a,b=linear_sum_assignment(-scores)
            similarities.append(float(scores[a,b].mean()))
        stability=min(similarities)
        if not (w.sum(axis=1)>0).all():
            raise ValueError('Active entities have degenerate NMF loadings')
        complete_w=np.zeros((len(matrix.entities),rank))
        complete_w[active]=w
        loadings=complete_w/np.maximum(complete_w.sum(axis=1,keepdims=True),1e-12)
        table(out/'loadings.tsv',[dict(entity=name,status='active' if active[i] else 'inactive',
              **{f'component_{j}':float(loadings[i,j]) for j in range(rank)})
              for i,name in enumerate(matrix.entities)],['entity','status',*[f'component_{j}' for j in range(rank)]])
        table(out/'components.tsv',[dict(block=block,feature=feature,**{f'component_{j}':float(h[j,i]) for j in range(rank)})
              for i,(block,feature) in enumerate(matrix.columns)],['block','feature',*[f'component_{j}' for j in range(rank)]])
        dump(out/'normalization.json',dict(features=scales,method='observed_95th_percentile_clipped',
              block_weights=settings.get('block_weights',{}),fitted_on_sha256=sha(features)))
        dump(out/'diagnostics.json',dict(**diagnostics,stability=stability,trials=[f[2] for f in fits],
             missing_entries=int((~matrix.observed).sum()),seed=seeds[best],
             inactive_entities=[name for i,name in enumerate(matrix.entities) if not active[i]],
             proposal_eligible=bool(stability>=float(settings.get('min_stability',.85)) and all(f[2]['converged'] for f in fits)),
             interpretation='Mixed-feature latent components are not taxa or a calibrated posterior; no automatic graph edits'))
        np.savez_compressed(out/'model.npz',w=complete_w,h=h,weights=weights,observed=matrix.observed,active=active)
