"""Run from repository root: snakemake -s workflows/reference.smk --cores 4."""
configfile: 'configs/workflow.yaml'

rule all:
    input:
        config['design'] + '/experiments.json'

rule reference:
    input:
        accessions=config['accessions'], settings='configs/reference.yaml'
    output:
        manifest=config['snapshot'] + '/manifest.json',
        fasta=config['snapshot'] + '/reference.fasta',
        metadata=config['snapshot'] + '/metadata.tsv'
    params:
        directory=config['snapshot']
    shell:
        'python -m nexvirome2 reference fetch --accessions {input.accessions:q} '
        '--config {input.settings:q} --output {params.directory:q}'

rule panel:
    input:
        fasta=rules.reference.output.fasta,
        metadata=rules.reference.output.metadata,
        manifest=rules.reference.output.manifest,
        settings='configs/panel.yaml'
    output:
        config['panel'] + '/panel.json'
    params:
        directory=config['panel']
    threads: config.get('threads', 4)
    shell:
        'python -m nexvirome2 panel build --reference {input.fasta:q} '
        '--metadata {input.metadata:q} --config {input.settings:q} --threads {threads} --output {params.directory:q}'

rule design:
    input:
        panel=rules.panel.output, settings=config['experiment_config']
    output:
        config['design'] + '/experiments.json'
    params:
        directory=config['design']
    shell:
        'python -m nexvirome2 prepare --panel {input.panel:q} '
        '--config {input.settings:q} --output {params.directory:q}'
