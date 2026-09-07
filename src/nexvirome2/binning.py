"""Conservative contig associations. Never concatenate sequences or infer taxonomy."""
import json
import itertools
from pathlib import Path
from collections import defaultdict
from .common import dump, rows, stage, table
from .coverage.math import cosine
from .coverage.matrix import read_matrix, aggregate


def build(settings,coverage_file,window_file,library_file,composition_file,output,profiles_file=None):
    import numpy as np
    inputs = [coverage_file,window_file,library_file,composition_file]+([profiles_file] if profiles_file else [])
    with stage(output,settings,inputs) as (out,manifest):
        regions,names,matrix = read_matrix(coverage_file,window_file)
        names,matrix = aggregate(matrix,names,library_file)
        comp = {r['window_id']:r['kmers'] for r in json.loads(Path(composition_file).read_text())}
        if set(comp)!={r['window_id'] for r in regions}:
            raise ValueError('Composition and coverage window IDs differ')
        by_contig = defaultdict(list)
        for i,r in enumerate(regions):
            by_contig[r['contig']].append(i)
        if len(by_contig)>settings.get('max_contigs',2000):
            raise ValueError('Candidate set too large for pairwise binning; select ambiguous neighborhoods')
        families = defaultdict(set)
        if profiles_file:
            for r in rows(profiles_file):
                if r.get('structure_family'):
                    families[r['contig']].add(r['structure_family'])
        features = {}
        for contig,indices in by_contig.items():
            profiles = matrix[indices]
            nonzero = profiles.sum(axis=1)>0
            if not nonzero.any():
                features[contig] = None
                continue
            p = profiles[nonzero].mean(axis=0)
            consistency = min(cosine(row,p) or 0 for row in profiles[nonzero])
            c = defaultdict(float)
            for i in indices:
                for word,value in comp[regions[i]['window_id']].items():
                    c[word] += value/len(indices)
            features[contig] = (p,c,consistency)
        parent = {n:n for n in by_contig}
        def root(n):
            while parent[n]!=n:
                parent[n]=parent[parent[n]]
                n=parent[n]
            return n
        edges = []
        for a,b in itertools.combinations(sorted(features),2):
            fa,fb = features[a],features[b]
            if fa is None or fb is None or min(fa[2],fb[2])<settings.get('min_internal_consistency',0.9):
                continue
            words = sorted(set(fa[1])|set(fb[1]))
            cov = cosine(fa[0],fb[0]) or 0
            composition = cosine([fa[1].get(w,0) for w in words],[fb[1].get(w,0) for w in words]) or 0
            union = families[a]|families[b]
            structure = len(families[a]&families[b])/len(union) if families[a] and families[b] else None
            weight = settings.get('structural_weight',0.5) if structure is not None else 0
            combined = (cov+composition+weight*(structure or 0))/(2+weight)
            if cov>=settings.get('min_coverage_similarity',0.95) and composition>=settings.get('min_composition_similarity',0.9) and combined>=settings.get('min_combined_similarity',0.9):
                edges.append({'a':a,'b':b,'coverage_similarity':cov,'composition_similarity':composition,'structural_jaccard':structure,'combined_similarity':combined})
                parent[root(b)] = root(a)
        groups = defaultdict(list)
        for name in sorted(parent):
            groups[root(name)].append(name)
        bins = []
        for i,members in enumerate(groups.values()):
            for name in members:
                bins.append({'contig':name,'bin':f'bin_{i:06d}',
                             'status':'associated' if len(members)>1 else 'unresolved' if features[name] is None or features[name][2]<settings.get('min_internal_consistency',0.9) else 'singleton'})
        table(out/'bins.tsv',bins,['contig','bin','status'])
        table(out/'edges.tsv',edges,['a','b','coverage_similarity','composition_similarity','structural_jaccard','combined_similarity'])
        manifest['interpretation'] = 'Connected candidate associations, not ordered genomes or demonstrated species boundaries'
