import warnings
from ..coverage.math import cosine

def factorize(matrix,settings):
    import numpy as np
    from sklearn.decomposition import NMF
    from scipy.optimize import linear_sum_assignment
    if matrix.ndim != 2 or not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError('Expected a finite nonnegative coverage matrix')
    if (matrix.sum(axis=0)>0).sum() < 3:
        raise ValueError('At least three biological samples must have nonzero coverage')
    active = matrix.sum(axis=1)>0
    x = matrix[active]
    if len(x)<2:
        raise ValueError('Need at least two nonzero windows')
    x = x/x.sum(axis=1,keepdims=True)
    if np.max(np.std(x,axis=0)) < 1e-8:
        raise ValueError('No variable abundance profiles; factorization is unidentifiable')
    seeds = settings.get('seeds',[17,29,43])
    if len(set(seeds)) < 2:
        raise ValueError('Use at least two distinct initialization seeds')
    max_rank = min(settings.get('max_rank',6),*x.shape)
    rank = settings.get('rank',2)
    ranks = range(2,max_rank+1) if rank == 'auto' else [int(rank)]
    trials = []
    for k in ranks:
        if not 2 <= k <= min(x.shape):
            raise ValueError('Rank outside matrix dimensions')
        fits = []
        for seed in seeds:
            model = NMF(n_components=k,init='random',random_state=seed,max_iter=settings.get('max_iter',2000),tol=1e-5)
            with warnings.catch_warnings(record=True) as caught:
                w = model.fit_transform(x)
            h = model.components_
            scale = h.sum(axis=1)
            h = h/np.maximum(scale[:,None],1e-12)
            w = w*scale
            error = float(np.linalg.norm(x-w@h)/np.linalg.norm(x))
            fits.append({'w':w,'h':h,'error':error,'seed':seed,'converged':model.n_iter_<model.max_iter,
                         'warnings':[str(c.message) for c in caught]})
        fits.sort(key=lambda f:f['error'])
        best = fits[0]
        stability = []
        for other in fits[1:]:
            similarity = np.array([[cosine(a,b) or 0 for b in other['h']] for a in best['h']])
            a,b = linear_sum_assignment(-similarity)
            stability.append(float(similarity[a,b].mean()))
        trials.append({'rank':k,'error':best['error'],'stability':min(stability),
                       'criterion':best['error']+settings.get('rank_penalty',0.02)*k,
                       'fits':fits,'best':best})
    winner = min(trials,key=lambda t:(t['criterion'],t['rank']))
    w = winner['best']['w']
    loadings = np.zeros((len(matrix),winner['rank']))
    loadings[active] = w/np.maximum(w.sum(axis=1,keepdims=True),1e-12)
    diagnostics = [{k:v for k,v in trial.items() if k not in ('fits','best')} |
                   {'runs':[{k:v for k,v in f.items() if k not in ('w','h')} for f in trial['fits']]} for trial in trials]
    return loadings,winner['best']['h'],winner,diagnostics,active


