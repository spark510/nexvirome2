import json
from collections import defaultdict
from ..common import dump,stage,table
from ..coverage.math import cosine
from ..coverage.matrix import read_matrix,aggregate
from .model import factorize

def fit(settings,coverage_file,window_file,library_file,output):
    import numpy as np
    with stage(output,settings,[coverage_file,window_file,library_file]) as (out,manifest):
        regions,names,matrix = read_matrix(coverage_file,window_file)
        biological,matrix = aggregate(matrix,names,library_file)
        loadings,profiles,winner,diagnostics,active = factorize(matrix,settings)
        fields = [f'factor_{i}' for i in range(loadings.shape[1])]
        table(out/'loadings.tsv',[{'window_id':r['window_id'],**dict(zip(fields,map(float,loading)))} for r,loading in zip(regions,loadings)],['window_id',*fields])
        table(out/'components.tsv',[{'component':fields[i],**dict(zip(biological,map(float,p)))} for i,p in enumerate(profiles)],['component',*biological])
        table(out/'profiles.tsv',[{'window_id':r['window_id'],**dict(zip(biological,map(float,p)))} for r,p in zip(regions,matrix)],['window_id',*biological])
        np.savez_compressed(out/'model.npz',components=profiles,loadings=loadings,samples=np.array(biological))
        dump(out/'diagnostics.json',{'selected_rank':winner['rank'],'stability':winner['stability'],'trials':diagnostics,
                                   'normalization':'normalized coverage averaged within biological sample, then row L1 for NMF'})
        by_contig = defaultdict(list)
        for i,r in enumerate(regions):
            by_contig[r['contig']].append(i)
        changes = []
        reliable = winner['stability'] >= settings.get('min_stability',0.85) and winner['best']['converged']
        for contig,indices in by_contig.items():
            indices.sort(key=lambda i:regions[i]['start'])
            for a,b in zip(indices,indices[1:]):
                if not active[a] or not active[b]:
                    continue
                left,right = regions[a],regions[b]
                # Do not bridge gaps introduced by a neighborhood-only window selection.
                if right['start'] > left['end']:
                    continue
                score = 1-(cosine(loadings[a],loadings[b]) or 0)
                cov_score = 1-(cosine(matrix[a],matrix[b]) or 0)
                high = min(loadings[a].max(),loadings[b].max()) >= settings.get('min_loading',0.7)
                if reliable and high and loadings[a].argmax()!=loadings[b].argmax() and score >= settings.get('min_change',0.5):
                    changes.append({'contig':contig,'start':min(left['start'],right['start']),
                                    'end':max(left['end'],right['end']),'method':'nmf','score':float(score),
                                    'detail':json.dumps({'left':left['window_id'],'right':right['window_id'],
                                                         'coverage_distance':cov_score,'stability':winner['stability']})})
        table(out/'proposals.tsv',changes,['contig','start','end','method','score','detail'])
        manifest.update(biological_samples=biological,rank=winner['rank'],reliable=reliable,
                        scientific_status='candidate evidence; not a chimera call')


