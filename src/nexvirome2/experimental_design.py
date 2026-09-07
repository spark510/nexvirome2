"""Biological abundance variation, separate from technical simulator seeds."""
from pathlib import Path
import math
import random
from .common import config, stage, dump, table, fasta


def multisample(settings, sources_file, output):
    spec=config(sources_file)
    sources=spec['sources']
    sample_count=int(settings.get('biological_samples',6))
    seeds=settings.get('technical_seeds',[17,29,43])
    base=float(settings.get('coverage',50))
    sigma=float(settings.get('log_abundance_sd',1))
    if len(sources)<2 or sample_count<3 or not seeds or len(set(seeds))!=len(seeds):
        raise ValueError('Need two sources, three biological samples and distinct technical seeds')
    if not math.isfinite(base) or base<=0 or not math.isfinite(sigma) or sigma<=0:
        raise ValueError('Coverage/abundance variation must be positive and finite')
    ids=[s['accession'] for s in sources]
    if len(set(ids))!=len(ids):
        raise ValueError('Duplicate source accession')
    paths=[(Path(sources_file).resolve().parent/s['fasta']).resolve() for s in sources]
    with stage(output,settings,[sources_file,*paths]) as (out,_):
        for source,path in zip(sources,paths):
            if source['accession'] not in fasta(path):
                raise ValueError('Source accession absent from source FASTA')
        rng=random.Random(int(settings.get('seed',17)))
        mode=settings.get('mode','independent')
        if mode not in ('independent','correlated','constant','sparse'):
            raise ValueError('Unknown abundance control mode')
        abundance=[]
        experiments={}
        for sample in range(sample_count):
            shared=rng.gauss(0,sigma)
            coverages=[]
            for i,source in enumerate(sources):
                noise=rng.gauss(0,sigma)
                log_value=0 if mode=='constant' else shared if mode=='correlated' else noise
                value=base*math.exp(log_value)
                if mode=='sparse' and (sample+i)%3==0:
                    value=0
                coverages.append(value)
                abundance.append(dict(biological_sample=f'b{sample}',accession=source['accession'],coverage=value))
            for seed in seeds:
                name=f'b{sample}_s{seed}'
                active=[{**source,'fasta':str(path),'coverage':value,'seed_offset':i}
                        for i,(source,path,value) in enumerate(zip(sources,paths,coverages)) if value>0]
                if not active:
                    raise ValueError('An empty biological sample requires a separately defined background control')
                dump(out/f'{name}.sources.json',{'sources':active})
                experiments[name]=dict(sources=str((out/f'{name}.sources.json').resolve()),seed=seed,
                     biological_sample=f'b{sample}',mode=mode)
        table(out/'abundance.tsv',abundance,['biological_sample','accession','coverage'])
        dump(out/'experiments.json',dict(experiments=experiments,biological_samples=sample_count,
             technical_seeds=seeds,mode=mode,constant_profile_control=mode in ('constant','correlated')))


def rank_paths(settings, features, output):
    """Simple weighted baseline on explicitly normalized [0,1] path evidence.

TSV: candidate,path,feature,value. Blank means missing, not a zero score.
No path ranking can authorize a sequence split or reconstruction by itself.
"""
    from collections import defaultdict
    from .common import rows
    weights=settings['weights']
    if not weights or any(not math.isfinite(float(v)) or float(v)<0 for v in weights.values()) or sum(map(float,weights.values()))<=0:
        raise ValueError('Feature weights must be finite, nonnegative, and not all zero')
    with stage(output,settings,[features]) as (out,_):
        groups=defaultdict(dict)
        for row in rows(features):
            key=(row['candidate'],row['path'])
            feature=row['feature']
            if not all(key) or feature not in weights or feature in groups[key]:
                raise ValueError('Duplicate/unknown path feature or empty ID')
            value=float(row['value']) if row['value'] else None
            if value is not None and not 0<=value<=1:
                raise ValueError('Path features must be normalized to [0,1]')
            groups[key][feature]=value
        result=[]
        for (candidate,path),values in groups.items():
            available={f:v for f,v in values.items() if v is not None and float(weights[f])>0}
            denominator=sum(float(weights[f]) for f in available)
            result.append(dict(candidate=candidate,path=path,
                 score=sum(float(weights[f])*v for f,v in available.items())/denominator if denominator else None,
                 observed_weight=denominator,missing=';'.join(sorted(set(weights)-set(available))),
                 decision='ranking_only'))
        result.sort(key=lambda r:(r['candidate'],r['score'] is None,-(r['score'] or 0),r['path']))
        table(out/'ranking.tsv',result,['candidate','path','score','observed_weight','missing','decision'])
