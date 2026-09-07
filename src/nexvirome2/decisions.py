"""Paired, independent-unit evidence gates; research runs never enable production."""
import hashlib
import json
import math
import random
from .common import rows, stage, dump, config, sha, table


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def freeze(settings, development, output):
    with stage(output,settings,[development]) as (out,_):
        groups = sorted({r['unit'] for r in rows(development)})
        if not groups or any(not x for x in groups):
            raise ValueError('Development units cannot be empty')
        parameters = settings['parameters']
        dump(out/'frozen.json',dict(schema=1,parameters=parameters,parameters_sha256=fingerprint(parameters),
             development_units=groups,development_sha256=sha(development)))


def paired_interval(values, seed=17, replicates=2000):
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError('Paired effects must be finite and nonempty')
    rng = random.Random(seed)
    effects = sorted(sum(rng.choices(values,k=len(values)))/len(values) for _ in range(replicates))
    return dict(mean=sum(values)/len(values),lower=effects[int(.025*replicates)],
                upper=effects[min(replicates-1,int(.975*replicates))])


def select(settings, evidence, frozen, output):
    """Evidence TSV: module,unit,baseline_error,error,baseline_recovery,recovery,resource_ratio.

Errors are rates per unit (not raw counts); recovery is a fraction. Unit IDs
must identify independent held-out reference/segment clusters or samples.
"""
    with stage(output,settings,[evidence,frozen]) as (out,_):
        lock = config(frozen)
        if lock.get('schema')!=1 or fingerprint(lock['parameters'])!=lock['parameters_sha256']:
            raise ValueError('Corrupt frozen parameters')
        if fingerprint(settings['parameters'])!=lock['parameters_sha256']:
            raise ValueError('Evaluation parameters differ from frozen parameters')
        records, groups = rows(evidence), {}
        for row in records:
            name,unit = row['module'],row['unit']
            if name not in lock['parameters']['modules'] or not unit:
                raise ValueError('Unknown module/empty evaluation unit')
            if unit in lock['development_units']:
                raise ValueError('Development/evaluation unit leakage')
            if unit in groups.setdefault(name,{}):
                raise ValueError('Duplicate module/unit evidence')
            parsed = {key:float(row[key]) for key in ('baseline_error','error','baseline_recovery','recovery','resource_ratio')}
            if not all(math.isfinite(v) for v in parsed.values()) or parsed['resource_ratio']<=0:
                raise ValueError('Invalid resource/effect value')
            if any(not 0<=parsed[k]<=1 for k in parsed if k!='resource_ratio'):
                raise ValueError('Error and recovery rates must be fractions')
            groups[name][unit] = parsed
        decisions = []
        parameters = lock['parameters']
        required = max(3,int(parameters.get('min_units',3)))
        for module,rule in parameters['modules'].items():
            data = list(groups.get(module,{}).values())
            reason, status, effects = 'no_heldout_evidence','research',{}
            if data:
                effects = {name:paired_interval([r[left]-r[right] for r in data])
                           for name,left,right in [('error_reduction','baseline_error','error'),
                                                   ('recovery_change','recovery','baseline_recovery')]}
                base = sum(r['baseline_error'] for r in data)/len(data)
                target = float(rule.get('min_error_reduction',0))
                relative = float(rule.get('min_relative_error_reduction',0))
                loss = float(rule.get('max_recovery_loss',.02))
                cost = float(rule.get('max_resource_ratio',2))
                if not all(math.isfinite(x) for x in (target,relative,loss,cost)) or not 0<=target<=1 or not 0<=relative<=1 or not 0<=loss<=1 or cost<=0:
                    raise ValueError('Invalid frozen adoption thresholds')
                if settings.get('evidence_kind')!='natural_heldout':
                    reason = 'fixtures_or_unverified_evidence'
                elif not isinstance(rule.get('implementation'),dict) or not rule['implementation'].get('command'):
                    reason = 'implementation_settings_not_frozen'
                elif len(data)<required:
                    reason = 'insufficient_independent_units'
                elif base==0:
                    reason = 'no_baseline_errors'
                elif effects['recovery_change']['upper'] < -loss or effects['error_reduction']['upper']<0:
                    status,reason = 'excluded','evidence_of_harm'
                elif max(r['resource_ratio'] for r in data)>cost:
                    status,reason = 'excluded','resource_limit_exceeded'
                elif (effects['error_reduction']['lower']>0 and effects['error_reduction']['lower']>=max(target,relative*base)
                      and effects['recovery_change']['lower']>=-loss):
                    status = 'conditional' if rule.get('requires') else 'default'
                    reason = 'passed_frozen_heldout_gate'
                else:
                    reason = 'benefit_or_noninferiority_not_established'
            decisions.append(dict(module=module,status=status,reason=reason,units=len(data),effects=effects,
                                  requires=rule.get('requires',[]),implementation_sha256=fingerprint(rule.get('implementation'))))
        dump(out/'decisions.json',dict(schema=1,decisions=decisions,parameters_sha256=lock['parameters_sha256'],
             evidence_sha256=sha(evidence),evidence_kind=settings.get('evidence_kind','unverified'),
             uncertainty='Paired bootstrap over independent units; small panels remain exploratory'))
        dump(out/'production.json',dict(schema=1,mandatory=['provenance','input_validation','uncertainty'],
             enabled=[d['module'] for d in decisions if d['status']=='default'],
             conditional={d['module']:d['requires'] for d in decisions if d['status']=='conditional'},
             research=[d['module'] for d in decisions if d['status']=='research'],
             excluded=[d['module'] for d in decisions if d['status']=='excluded'],
             implementations={name:rule.get('implementation') for name,rule in parameters['modules'].items()},
             decision_sha256=sha(out/'decisions.json')))
        table(out/'decisions.tsv',decisions,['module','status','reason','units'])
