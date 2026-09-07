"""Select whole reference-mask policies on paired development experiments only."""
from collections import defaultdict
from pathlib import Path
import math
from .common import config, rows, stage, dump, sha, table
from .masking import MaskSet
from .decisions import paired_interval, fingerprint


def optimize(settings, reference, candidates_file, measurements, output):
    """Candidates JSON: conditions mapping and policies mapping (unmasked:null).

    Measurements TSV: condition,policy,unit,queries_sha256,wrong,correct,
    unclassified,retained_breadth. Each policy must have identical units and query
    fingerprints within its condition. Breadth must be measured externally on the
    unmasked coordinate system, not substituted by assigned-read fraction.
    """
    spec=config(candidates_file)
    if settings.get('partition')!='development':
        raise ValueError('Policy optimization is restricted to development data')
    if spec['policies'].get('unmasked','missing') is not None:
        raise ValueError('An explicit unmasked baseline is required')
    paths={name:(Path(candidates_file).resolve().parent/path).resolve()
           for name,path in spec['policies'].items() if path is not None}
    if len(paths)+1!=len(spec['policies']):
        raise ValueError('Only the unmasked baseline may have a null policy')
    constraints={name:float(settings.get(name,.02)) for name in ('max_correct_loss','max_breadth_loss')}
    if any(not math.isfinite(v) or not 0<=v<=1 for v in constraints.values()):
        raise ValueError('Policy loss constraints must be fractions')
    with stage(output,settings,[reference,candidates_file,measurements,*paths.values()]) as (out,_):
        masks={name:MaskSet.load(path,reference) for name,path in paths.items()}
        for condition in spec['conditions'].values():
            if not isinstance(condition,dict) or not all(k in condition for k in ('rank','read_length','library','database_sha256')):
                raise ValueError('Conditions require rank, read_length, library and database_sha256')
            if int(condition['read_length'])<1 or condition['database_sha256']!=sha(reference):
                raise ValueError('Condition reference checksum/read length mismatch')
        groups=defaultdict(lambda:defaultdict(dict))
        for row in rows(measurements):
            condition,policy,unit=row['condition'],row['policy'],row['unit']
            if condition not in spec['conditions'] or policy not in spec['policies'] or not unit:
                raise ValueError('Unknown condition/policy or empty independent unit')
            if unit in groups[condition][policy]:
                raise ValueError('Duplicate development condition/policy/unit')
            digest=row['queries_sha256']
            if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
                raise ValueError('Invalid query fingerprint')
            parsed={key:float(row[key]) for key in ('wrong','correct','unclassified','retained_breadth')}
            if any(not math.isfinite(v) or not 0<=v<=1 for v in parsed.values()):
                raise ValueError('Empirical measurements must be finite fractions')
            if sum(parsed[k] for k in ('wrong','correct','unclassified'))>1+1e-9:
                raise ValueError('Classification fractions exceed the all-query denominator')
            groups[condition][policy][unit]={**parsed,'queries_sha256':digest}
        if set(groups)!=set(spec['conditions']):
            raise ValueError('Every configured condition needs measurements')
        audit,selection=[],{}
        for condition,policies in sorted(groups.items()):
            applicable={'unmasked',*[name for name,mask in masks.items()
                                     if mask.policy['rank']==spec['conditions'][condition]['rank']]}
            if set(policies)!=applicable:
                raise ValueError('Every condition must evaluate every candidate policy at its rank')
            baseline=policies['unmasked']
            if len(baseline)<max(3,int(settings.get('min_units',3))):
                raise ValueError('At least three independent development units are required')
            eligible=[]
            for name,values in sorted(policies.items()):
                if values.keys()!=baseline.keys() or any(values[u]['queries_sha256']!=baseline[u]['queries_sha256'] for u in baseline):
                    raise ValueError('Policies must use identical paired units and queries')
                if name=='unmasked':
                    continue
                if masks[name].policy['rank']!=spec['conditions'][condition]['rank']:
                    raise ValueError('Candidate mask rank differs from experimental condition')
                effects={key:paired_interval([values[u][key]-baseline[u][key] for u in sorted(baseline)])
                         for key in ('wrong','correct','unclassified','retained_breadth')}
                passed=(effects['wrong']['upper']<0 and
                        effects['correct']['lower']>=-constraints['max_correct_loss'] and
                        effects['retained_breadth']['lower']>=-constraints['max_breadth_loss'])
                audit.append(dict(condition=condition,policy=name,eligible=passed,effects=effects))
                if passed:
                    eligible.append((effects['wrong']['mean'],-effects['correct']['mean'],name))
            chosen=min(eligible)[2] if eligible else 'unmasked'
            condition_hash=fingerprint(spec['conditions'][condition])
            saved=None
            if chosen!='unmasked':
                saved=f'{condition_hash}.mask.json'
                old=masks[chosen]
                policy={**old.policy,'conditions':spec['conditions'][condition],
                        'selection_partition':'development','selection_measurements_sha256':sha(measurements)}
                MaskSet.create(reference,old.intervals,policy).save(out/saved)
            selection[condition]=dict(policy=chosen,mask=saved,conditions=spec['conditions'][condition],
                                       development_units=sorted(baseline),condition_sha256=condition_hash)
        dump(out/'selection.json',dict(schema=1,partition='development',selection=selection,
             measurements_sha256=sha(measurements),production_approved=False,
             note='Development policy selection is subject to independent heldout evaluation'))
        dump(out/'candidate_effects.json',audit)
        table(out/'selection.tsv',[dict(condition=c,policy=v['policy'],mask=v['mask']) for c,v in selection.items()],
              ['condition','policy','mask'])
        table(out/'development_units.tsv',[dict(unit=u) for u in sorted({u for v in selection.values() for u in v['development_units']})],['unit'])
