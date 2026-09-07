"""CLI registration for reference masking and evidence-gated adoption."""


COMMANDS = {
    'taxonomy': {'build':['mapping','nodes']},
    'sharing': {'build':['reference','taxonomy'], 'align':['reference','taxonomy']},
    'mask': {'import':['reference','source'], 'export':['reference','mask'],
             'propose':['reference','atlas'],
             'optimize':['reference','candidates','measurements'],
             'evaluate':['queries','reference','taxonomy','truth','arms']},
    'classification': {'run':['queries','reference','taxonomy'],
                       'evaluate':['predictions','truth','metadata']},
    'pipeline': {'freeze':['development'], 'select':['evidence','frozen'],
                 'execute':['production','decisions','recipes']},
    'evidence': {'link':['reference','mask','annotations'], 'domains':['reference','orfs','hits'],
                 'joint':['reference','observations']},
    'discovery': {'summarize':['profiles']},
    'catalog': {'build':['inputs'], 'split':['membership'], 'cluster':['contigs']},
    'design': {'multisample':['sources']},
    'paths': {'rank':['features']},
    'quality': {'genomad':['contigs','summary']},
    'multiview': {'fit':['features'],'propose':['model','coordinates']},
    'reads': {'prepare':['r1','r2']},
    'conserved': {'reuse':['masked','unmasked']},
    'coassembly': {'prepare':['samples','membership']},
}


def register(commands):
    for group,actions in COMMANDS.items():
        sub = commands.add_parser(group).add_subparsers(dest='action',required=True)
        for action,fields in actions.items():
            parser = sub.add_parser(action)
            for field in ['config','output',*fields]:
                parser.add_argument('--'+field,required=True)
            if group=='taxonomy':
                parser.add_argument('--merged')
                parser.add_argument('--deleted')
            if group=='classification' and action=='run':
                parser.add_argument('--mask')
                parser.add_argument('--hits')
            if group=='discovery':
                parser.add_argument('--motifs')
                parser.add_argument('--context')
            if group=='sharing':
                parser.add_argument('--metadata' if action=='build' else '--alignments')
            if group=='reads':
                parser.add_argument('--host-calls')
                parser.add_argument('--host-metadata')
            if group=='catalog' and action=='cluster':
                parser.add_argument('--alignments')


def dispatch(args,settings):
    if args.command=='taxonomy':
        from .taxonomy import build
        build(settings,args.mapping,args.nodes,args.output,args.merged,args.deleted)
    elif args.command=='sharing':
        if args.action=='build':
            from .sharing import build
            build(settings,args.reference,args.taxonomy,args.output,args.metadata)
        else:
            from .sharing_alignment import build
            build(settings,args.reference,args.taxonomy,args.output,args.alignments)
    elif args.command=='mask':
        from .masking import import_mask,export_mask,propose
        from .mask_benchmark import run
        if args.action=='import':
            import_mask(settings,args.reference,args.source,args.output)
        elif args.action=='export':
            export_mask(settings,args.reference,args.mask,args.output)
        elif args.action=='propose':
            propose(settings,args.reference,args.atlas,args.output)
        elif args.action=='optimize':
            from .mask_optimization import optimize
            optimize(settings,args.reference,args.candidates,args.measurements,args.output)
        else:
            run(settings,args.queries,args.reference,args.taxonomy,args.truth,args.arms,args.output)
    elif args.command=='classification':
        from .classification import classify,evaluate
        if args.action=='run':
            classify(settings,args.queries,args.reference,args.taxonomy,args.output,args.mask,args.hits)
        else:
            evaluate(settings,args.predictions,args.truth,args.output,args.metadata)
    elif args.command=='pipeline':
        from .decisions import freeze,select
        if args.action=='freeze':
            freeze(settings,args.development,args.output)
        elif args.action=='select':
            select(settings,args.evidence,args.frozen,args.output)
        else:
            from .production import execute
            execute(settings,args.production,args.decisions,args.recipes,args.output)
    elif args.command=='evidence':
        from .integrated_evidence import link,domains,joint
        if args.action=='link':
            link(settings,args.reference,args.mask,args.annotations,args.output)
        elif args.action=='domains':
            domains(settings,args.reference,args.orfs,args.hits,args.output)
        else:
            joint(settings,args.reference,args.observations,args.output)
    elif args.command=='discovery':
        from .integrated_evidence import discovery
        discovery(settings,args.profiles,args.output,args.motifs,args.context)
    elif args.command=='catalog':
        from .catalog import build,split
        if args.action=='build':
            build(settings,args.inputs,args.output)
        elif args.action=='split':
            split(settings,args.membership,args.output)
        else:
            from .catalog_clustering import cluster
            cluster(settings,args.contigs,args.output,args.alignments)
    elif args.command=='design':
        from .experimental_design import multisample
        multisample(settings,args.sources,args.output)
    elif args.command=='paths':
        from .experimental_design import rank_paths
        rank_paths(settings,args.features,args.output)
    elif args.command=='quality':
        from .quality_results import genomad
        genomad(settings,args.contigs,args.summary,args.output)
    elif args.command=='multiview':
        from .multiview import fit,propose
        if args.action=='fit':
            fit(settings,args.features,args.output)
        else:
            propose(settings,args.model,args.coordinates,args.output)
    elif args.command=='reads':
        from .read_preparation import prepare
        prepare(settings,args.r1,args.r2,args.output,args.host_calls,args.host_metadata)
    elif args.command=='conserved':
        from .conserved_reuse import link
        link(settings,args.masked,args.unmasked,args.output)
    elif args.command=='coassembly':
        from .coassembly import prepare
        prepare(settings,args.samples,args.membership,args.output)
