"""Optional common-target multi-sample analysis and correction ablations."""
from pathlib import Path
from nexvirome2.common import config as read_config, rows, dump, table, stage
from nexvirome2.searches import artifact_files

configfile: 'configs/extensions_workflow.yaml'
ROOT = config['output']
VARIANTS = ['baseline','simple','nmf','structure','combined']
if config.get('multiview_features'):
    VARIANTS.append('multiview')
SAMPLE_INPUTS = [str((Path(config['samples']).resolve().parent/r[k]).resolve())
                 for r in rows(config['samples']) for k in (['sam'] if r.get('sam') else ['r1','r2'])]
DB_INPUTS = [str(p) for key in ['sequence_database','structure_database','prostt5_model'] for p in artifact_files(config[key])]

rule all:
    input:
        ROOT+'/summary/ablations.tsv',
        [ROOT+'/accuracy/metrics.tsv'] if config.get('truth_sources') else []

rule neighborhoods:
    input: contigs=config['contigs'], graph=config['graph'], paths=config['paths'], settings='configs/features.yaml'
    output: ROOT+'/neighborhoods/regions.tsv'
    run:
        from nexvirome2.features import neighborhoods
        neighborhoods(read_config(input.settings),input.contigs,input.graph,input.paths,str(Path(output[0]).parent))

rule coverage:
    input:
        contigs=config['contigs'], samples=config['samples'], reads=SAMPLE_INPUTS,
        regions=ROOT+'/neighborhoods/regions.tsv' if config.get('local_only',False) else [],
        settings='configs/features.yaml'
    output:
        matrix=ROOT+'/features/normalized_coverage.tsv', windows=ROOT+'/features/windows.tsv',
        libraries=ROOT+'/features/libraries.tsv', composition=ROOT+'/features/composition.json'
    threads: config.get('threads',4)
    run:
        from nexvirome2.features import coverage
        coverage(read_config(input.settings)|{'threads':threads},input.contigs,input.samples,ROOT+'/features',
                 str(input.regions) if config.get('local_only',False) else None)

rule simple:
    input:
        matrix=rules.coverage.output.matrix, windows=rules.coverage.output.windows,
        libraries=rules.coverage.output.libraries, composition=rules.coverage.output.composition,
        settings='configs/features.yaml'
    output: ROOT+'/simple/proposals.tsv'
    run:
        from nexvirome2.features import propose
        propose(read_config(input.settings),input.matrix,input.windows,input.libraries,input.composition,ROOT+'/simple')

rule nmf:
    input:
        matrix=rules.coverage.output.matrix, windows=rules.coverage.output.windows,
        libraries=rules.coverage.output.libraries, settings='configs/latent.yaml'
    output: proposals=ROOT+'/nmf/proposals.tsv', diagnostics=ROOT+'/nmf/diagnostics.json'
    run:
        from nexvirome2.latent import fit
        fit(read_config(input.settings),input.matrix,input.windows,input.libraries,ROOT+'/nmf')

rule orfs:
    input:
        contigs=config['contigs'], settings='configs/structure.yaml', imported=config.get('orf_coordinates',[]),
        regions=ROOT+'/neighborhoods/regions.tsv' if config.get('local_only',False) else []
    output: proteins=ROOT+'/orfs/proteins.faa', coordinates=ROOT+'/orfs/orfs.tsv'
    run:
        from nexvirome2.orfs import predict
        predict(read_config(input.settings),input.contigs,ROOT+'/orfs',config.get('orf_coordinates'),
                str(input.regions) if config.get('local_only',False) else None)

rule search:
    input: proteins=rules.orfs.output.proteins, artifacts=DB_INPUTS, settings='configs/structure.yaml'
    output: ROOT+'/search/{backend}/hits.tsv'
    wildcard_constraints: backend='sequence|structure'
    threads: config.get('threads',4)
    run:
        from nexvirome2.application import run_search, SearchRequest
        structural = wildcards.backend=='structure'
        request = SearchRequest(input.proteins,config['structure_database' if structural else 'sequence_database'],
                                'foldseek' if structural else 'mmseqs',model=config['prostt5_model'] if structural else None,
                                cache=str(Path(config.get('cache',ROOT+'/cache'))/(wildcards.backend+'.sqlite')))
        run_search(read_config(input.settings)|{'threads':threads},request,str(Path(output[0]).parent))

rule annotate:
    input:
        coordinates=rules.orfs.output.coordinates, sequence=ROOT+'/search/sequence/hits.tsv',
        structure=ROOT+'/search/structure/hits.tsv', metadata=config['target_metadata'], settings='configs/structure.yaml'
    output: proposals=ROOT+'/annotation/proposals.tsv', profiles=ROOT+'/annotation/orf_profiles.tsv', contigs=ROOT+'/annotation/contigs.tsv'
    run:
        from nexvirome2.structural import annotation
        annotation(read_config(input.settings),input.coordinates,input.metadata,ROOT+'/annotation',input.sequence,input.structure)

rule bins:
    input:
        matrix=rules.coverage.output.matrix, windows=rules.coverage.output.windows,
        libraries=rules.coverage.output.libraries, composition=rules.coverage.output.composition,
        profiles=rules.annotate.output.profiles, settings='configs/binning.yaml'
    output: ROOT+'/bins/bins.tsv'
    run:
        from nexvirome2.binning import build
        build(read_config(input.settings),input.matrix,input.windows,input.libraries,input.composition,ROOT+'/bins',input.profiles)

def proposals_for(w):
    return {'baseline':[], 'simple':[ROOT+'/simple/proposals.tsv'], 'nmf':[ROOT+'/nmf/proposals.tsv'],
            'structure':[ROOT+'/annotation/proposals.tsv'],
            'multiview':[ROOT+'/multiview/proposals/proposals.tsv'],
            'combined':[ROOT+'/nmf/proposals.tsv',ROOT+'/annotation/proposals.tsv']}[w.variant]

rule multiview_fit:
    input:
        features=lambda w: config['multiview_features'], settings='configs/multiview.yaml'
    output:
        loadings=ROOT+'/multiview/fit/loadings.tsv', diagnostics=ROOT+'/multiview/fit/diagnostics.json'
    run:
        from nexvirome2.multiview import fit
        fit(read_config(input.settings),input.features,ROOT+'/multiview/fit')

rule multiview_proposals:
    input:
        loadings=rules.multiview_fit.output.loadings, diagnostics=rules.multiview_fit.output.diagnostics,
        coordinates=lambda w: config['multiview_coordinates'], settings='configs/multiview.yaml'
    output: ROOT+'/multiview/proposals/proposals.tsv'
    run:
        from nexvirome2.multiview import propose
        propose(read_config(input.settings),ROOT+'/multiview/fit',input.coordinates,ROOT+'/multiview/proposals')

rule correct:
    input:
        contigs=config['contigs'], graph=config['graph'], paths=config['paths'],
        r1=config['r1'], r2=config['r2'], proposals=proposals_for, settings='configs/correction.yaml'
    output:
        summary=ROOT+'/correction/{variant}/summary.json',
        corrected=ROOT+'/correction/{variant}/corrected.fasta',
        random=ROOT+'/correction/{variant}/random_control.fasta'
    wildcard_constraints: variant='baseline|simple|nmf|structure|combined|multiview'
    threads: config.get('threads',4)
    run:
        from nexvirome2.application import run_correction, CorrectionRequest, ablation_settings
        settings = ablation_settings(read_config(input.settings),wildcards.variant,threads)
        request = CorrectionRequest(input.contigs,input.graph,input.paths,input.r1,input.r2,tuple(input.proposals))
        run_correction(settings,request,str(Path(output.summary).parent))

def accuracy_contigs(w):
    if w.arm=='original':
        return config['contigs']
    variant,kind=w.arm.rsplit('_',1)
    return ROOT+'/correction/'+variant+('/corrected.fasta' if kind=='corrected' else '/random_control.fasta')

ACCURACY_ARMS=['original']+[variant+'_'+kind for variant in VARIANTS for kind in ['corrected','random']]

rule evaluate_ablation:
    input:
        contigs=accuracy_contigs,
        sources=lambda w: config['truth_sources'],
        settings='configs/evaluation.yaml'
    output: ROOT+'/evaluation/{arm}/metrics.json'
    wildcard_constraints: arm='original|(?:baseline|simple|nmf|structure|combined|multiview)_(?:corrected|random)'
    threads: config.get('threads',4)
    run:
        from nexvirome2.evaluation import evaluate
        evaluate(read_config(input.settings)|{'threads':threads},input.contigs,input.sources,str(Path(output[0]).parent))

rule accuracy_summary:
    input: expand(ROOT+'/evaluation/{arm}/metrics.json',arm=ACCURACY_ARMS)
    output: ROOT+'/accuracy/metrics.tsv'
    run:
        with stage(ROOT+'/accuracy',dict(config),list(input)) as (out,manifest):
            metrics=[]
            source_metrics=[]
            for arm,path in zip(ACCURACY_ARMS,input):
                for cutoff,values in read_config(path).items():
                    metrics.append(dict(arm=arm,min_length=cutoff,**{k:v for k,v in values.items() if k!='per_source'}))
                    source_metrics.extend(dict(arm=arm,min_length=cutoff,accession=acc,**value)
                                          for acc,value in values.get('per_source',{}).items())
            table(out/'metrics.tsv',metrics,['arm','min_length','contigs','chimeric_contigs','wrong_connections',
                  'chimera_rate','genome_recovery','correctly_reconstructed_fraction','ambiguous_contigs'])
            table(out/'sources.tsv',source_metrics,['arm','min_length','accession','recovered_bp','recovery',
                  'longest_aligned_block_bp','longest_consistent_block_bp'])
            dump(out/'interpretation.json',{'arms':ACCURACY_ARMS,'truth_used_only_for_evaluation':True,
                 'random_controls':'Each variant has its own matched split-count control',
                 'generalization_claim':False})

rule summary:
    input:
        corrections=expand(ROOT+'/correction/{variant}/summary.json',variant=VARIANTS),
        bins=rules.bins.output[0], annotations=rules.annotate.output.contigs, diagnostics=rules.nmf.output.diagnostics,
        simple=rules.simple.output[0], nmf=rules.nmf.output.proposals, structure=rules.annotate.output.proposals
    output: ROOT+'/summary/ablations.tsv'
    run:
        import json
        with stage(ROOT+'/summary',dict(config),list(input)) as (out,manifest):
            result = [dict(variant=variant,**json.loads(Path(path).read_text())) for variant,path in zip(VARIANTS,input.corrections)]
            table(out/'ablations.tsv',result,['variant','split_count','unresolved','truth_used'])
            dump(out/'interpretation.json',{'status':'Counts only; use independent truth-based evaluate for accuracy and recovery',
                'proposals':{k:len(rows(getattr(input,k))) for k in ['simple','nmf','structure']},
                'combined_rule':'read-supported split AND nmf AND structural_lineage'})
