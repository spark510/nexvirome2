"""Reference sharing -> policy and matched control -> identical-query evaluation."""
from pathlib import Path
import json

configfile: 'configs/masking_workflow.yaml'
ROOT = config['output']

rule all:
    input: ROOT + '/benchmark/comparison.json'

rule sharing_atlas:
    input:
        reference=config['reference'], taxonomy=config['taxonomy'], settings=config['settings'],
        metadata=[config['metadata']] if config.get('metadata') else []
    output: ROOT + '/sharing/atlas.tsv'
    params: directory=ROOT + '/sharing'
    run:
        from nexvirome2.common import config as read_config
        from nexvirome2.sharing import build
        build(read_config(input.settings),input.reference,input.taxonomy,params.directory,
              str(input.metadata[0]) if input.metadata else None)

rule mask_policy:
    input:
        reference=config['reference'], atlas=rules.sharing_atlas.output, settings=config['settings']
    output:
        proposed=ROOT + '/policy/mask.json', random=ROOT + '/policy/random_mask.json'
    params: directory=ROOT + '/policy'
    shell:
        'python -m nexvirome2 mask propose --reference {input.reference:q} --atlas {input.atlas:q} '
        '--config {input.settings:q} --output {params.directory:q}'

rule mask_arms:
    input:
        proposed=rules.mask_policy.output.proposed, random=rules.mask_policy.output.random,
        legacy=[config['legacy_mask']] if config.get('legacy_mask') else []
    output: ROOT + '/arms.json'
    run:
        arms = {'unmasked':None,'shared':str(Path(input.proposed).resolve()),
                'random':str(Path(input.random).resolve())}
        if input.legacy:
            arms['legacy'] = str(Path(input.legacy[0]).resolve())
        Path(output[0]).write_text(json.dumps(arms,indent=2)+'\n',encoding='utf-8')

rule mask_benchmark:
    input:
        queries=config['queries'], truth=config['truth'], reference=config['reference'],
        taxonomy=config['taxonomy'], arms=rules.mask_arms.output, settings=config['settings'],
        proposed=rules.mask_policy.output.proposed, random=rules.mask_policy.output.random,
        legacy=[config['legacy_mask']] if config.get('legacy_mask') else []
    output: ROOT + '/benchmark/comparison.json'
    params: directory=ROOT + '/benchmark'
    shell:
        'python -m nexvirome2 mask evaluate --queries {input.queries:q} --reference {input.reference:q} '
        '--taxonomy {input.taxonomy:q} --truth {input.truth:q} --arms {input.arms:q} '
        '--config {input.settings:q} --output {params.directory:q}'
