"""Execute only evidence-approved modules using frozen implementations."""
from pathlib import Path
import re
from .common import config,stage,dump,sha
from .decisions import fingerprint


ALLOWED = {('multiview','fit'),('multiview','propose'),('assemble',),('correct',),('features','coverage'),('features','neighborhoods'),
           ('features','propose'),('latent','fit'),('orfs','predict'),('search','sequence'),
           ('search','structure'),('structure','annotate'),('bins','build'),('qc','run'),
           ('mask','export'),('classification','run'),('catalog','build'),('discovery','summarize')}

DATA_ARGUMENTS = {
    ('multiview','fit'): {'features'}, ('multiview','propose'): {'model','coordinates'},
    ('assemble',): {'assembler','r1','r2'},
    ('correct',): {'contigs','graph','paths','r1','r2','proposals'},
    ('features','coverage'): {'contigs','samples','regions'},
    ('features','neighborhoods'): {'contigs','graph','paths'},
    ('features','propose'): {'coverage','windows','libraries','composition'},
    ('latent','fit'): {'coverage','windows','libraries'},
    ('orfs','predict'): {'contigs','imported','regions'},
    ('search','sequence'): {'proteins','database','cache'},
    ('search','structure'): {'proteins','database','model','structure-map','cache'},
    ('structure','annotate'): {'orfs','metadata','sequence-hits','structure-hits'},
    ('bins','build'): {'coverage','windows','libraries','composition','profiles'},
    ('qc','run'): {'tool','contigs','database'}, ('mask','export'): {'reference','mask'},
    ('classification','run'): {'queries','reference','taxonomy','mask','hits'},
    ('catalog','build'): {'inputs'}, ('discovery','summarize'): {'profiles','motifs','context'},
}


def execute(settings, production_file, decisions_file, recipes_file, output):
    """Recipes list {module,arguments,depends_on}; configuration is frozen in policy.

`${module}/artifact` resolves a completed predecessor. Conditions are explicit
deployment facts in settings.conditions; absent conditions disable that module.
Simulation/evaluation/truth operations are not part of this execution interface.
"""
    with stage(output,settings,[production_file,decisions_file,recipes_file]) as (out,_):
        policy,decisions,recipes=config(production_file),config(decisions_file),config(recipes_file)
        if policy.get('schema')!=1 or policy['decision_sha256']!=sha(decisions_file):
            raise ValueError('Production decision checksum mismatch')
        expected_enabled={r['module'] for r in decisions['decisions'] if r['status']=='default'}
        expected_conditional={r['module']:r['requires'] for r in decisions['decisions'] if r['status']=='conditional'}
        implementation_hashes={r['module']:r['implementation_sha256'] for r in decisions['decisions']}
        if set(policy['enabled'])!=expected_enabled or policy['conditional']!=expected_conditional:
            raise ValueError('Production activation differs from evidence decision')
        names=[r['module'] for r in recipes['stages']]
        if len(set(names))!=len(names) or any(not re.fullmatch('[a-z0-9_-]+',n) for n in names):
            raise ValueError('Duplicate/unsafe recipe module name')
        active=set(policy['enabled'])|set(policy['conditional'])
        if active-set(names):
            raise ValueError('An approved module has no execution recipe')
        completed,logs={},[]
        facts=settings.get('conditions',{})
        for recipe in recipes['stages']:
            name=recipe['module']
            if name not in active:
                logs.append(dict(module=name,status='skipped',reason='not_approved'))
                continue
            requirements=policy['conditional'].get(name,[])
            missing=[condition for condition in requirements if facts.get(condition) is not True]
            if missing:
                logs.append(dict(module=name,status='skipped',reason='unmet_conditions',conditions=missing))
                continue
            if any(dep not in completed for dep in recipe.get('depends_on',[])):
                raise ValueError('Recipe dependencies must precede module and complete successfully')
            implementation=policy['implementations'][name]
            if fingerprint(implementation)!=implementation_hashes[name]:
                raise ValueError('Implementation differs from frozen evaluation decision')
            if not isinstance(implementation,dict) or tuple(implementation['command']) not in ALLOWED:
                raise ValueError('Unapproved production operation')
            command=implementation['command']
            module_settings=implementation.get('settings',{})
            def resolve(value):
                text=str(value)
                match=re.match(r'^\$\{([a-z0-9_-]+)\}/(.+)$',text)
                if match:
                    if match[1] not in completed:
                        raise ValueError('Referenced module has not completed')
                    base=completed[match[1]].resolve()
                    result=(base/match[2]).resolve()
                    if not result.is_relative_to(base):
                        raise ValueError('Artifact path escapes predecessor output')
                    return str(result)
                return text
            config_file=out/f'{name}.config.json'
            dump(config_file,module_settings)
            argv=[*command,'--config',str(config_file),'--output',str(out/name)]
            for key,value in recipe.get('arguments',{}).items():
                if key not in DATA_ARGUMENTS[tuple(command)]:
                    raise ValueError('Reserved or invalid production argument')
                argv.append('--'+key)
                values=[resolve(v) for v in (value if isinstance(value,list) else [value])]
                if any(v.startswith('--') for v in values):
                    raise ValueError('Production input values must not contain option tokens')
                argv.extend(values)
            from .cli import main
            try:
                result=main(argv)
            except SystemExit as exc:
                raise ValueError(f'Invalid recipe arguments: {name}') from exc
            if result:
                raise RuntimeError(f'Production module failed: {name}')
            module_manifest=config(out/name/'manifest.json')
            if module_manifest.get('status')!='complete':
                raise RuntimeError('Production module did not complete')
            if fingerprint(config(out/name/'config.json'))!=fingerprint(module_settings):
                raise RuntimeError('Effective module settings differ from frozen implementation')
            completed[name]=out/name
            logs.append(dict(module=name,status='complete',implementation_sha256=fingerprint(implementation)))
            dump(out/'execution.json',dict(stages=logs))
        dump(out/'execution.json',dict(stages=logs,truth_used=False))
