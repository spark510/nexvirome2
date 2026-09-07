"""Compare explicitly selected pools with matched whole per-sample assemblies."""
from pathlib import Path
from nexvirome2.common import rows, config as read_config, stage, table, fasta

ROOT=config['output']
SAMPLES={row['sample']:row for row in rows(config['samples'])}
GROUPS=sorted({row['group'] for row in rows(config['membership'])})
ASSEMBLERS=config.get('assemblers',['metaspades','megahit'])
if not set(ASSEMBLERS)<= {'metaspades','megahit'} or not ASSEMBLERS:
    raise ValueError('Unsupported/empty coassembly assembler selection')
BASE=Path(config['samples']).resolve().parent

rule all:
    input: ROOT+'/comparison/assemblies.tsv'

rule pool_reads:
    input:
        samples=config['samples'], membership=config['membership'],
        reads=[str((BASE/row[k]).resolve()) for row in SAMPLES.values() for k in ('r1','r2')]
    output:
        manifest=ROOT+'/pools/manifest.json',
        reads=expand(ROOT+'/pools/{group}/reads_R{mate}.fastq.gz',group=GROUPS,mate=[1,2])
    run:
        from nexvirome2.coassembly import prepare
        prepare(dict(config),input.samples,input.membership,ROOT+'/pools')

rule pooled_assembly:
    input:
        r1=ROOT+'/pools/{group}/reads_R1.fastq.gz',r2=ROOT+'/pools/{group}/reads_R2.fastq.gz',
        manifest=rules.pool_reads.output.manifest, settings='configs/assembly.yaml'
    output: ROOT+'/pooled/{group}/{assembler}/result.json'
    wildcard_constraints: assembler='metaspades|megahit'
    threads: config.get('threads',4)
    run:
        from nexvirome2.assembly import assemble
        assemble(read_config(input.settings)|{'threads':threads},wildcards.assembler,input.r1,input.r2,str(Path(output[0]).parent))

rule per_sample_assembly:
    input:
        r1=lambda w: str((BASE/SAMPLES[w.sample]['r1']).resolve()),
        r2=lambda w: str((BASE/SAMPLES[w.sample]['r2']).resolve()), settings='configs/assembly.yaml'
    output: ROOT+'/per_sample/{sample}/{assembler}/result.json'
    wildcard_constraints: assembler='metaspades|megahit'
    threads: config.get('threads',4)
    run:
        from nexvirome2.assembly import assemble
        assemble(read_config(input.settings)|{'threads':threads},wildcards.assembler,input.r1,input.r2,str(Path(output[0]).parent))

rule comparison:
    input:
        expand(ROOT+'/pooled/{group}/{assembler}/result.json',group=GROUPS,assembler=ASSEMBLERS),
        expand(ROOT+'/per_sample/{sample}/{assembler}/result.json',sample=sorted(SAMPLES),assembler=ASSEMBLERS)
    output: ROOT+'/comparison/assemblies.tsv'
    run:
        artifacts=[Path(path).parent/read_config(path)['contigs'] for path in input]
        with stage(ROOT+'/comparison',dict(config),[*input,*artifacts]) as (out,manifest):
            result=[]
            for path,contigs in zip(input,artifacts):
                path=Path(path)
                lengths=[len(seq) for seq in fasta(contigs).values()]
                result.append(dict(mode=path.parents[2].name,sample_or_group=path.parents[1].name,
                    assembler=path.parent.name,contigs=len(lengths),bases=sum(lengths),longest=max(lengths,default=0)))
            table(out/'assemblies.tsv',result,['mode','sample_or_group','assembler','contigs','bases','longest'])
            manifest['interpretation']='Assembly size comparison only; source-aware accuracy and read-budget matching require independent evaluation'
