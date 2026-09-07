"""Explicit CLI: correction does not accept simulation truth."""
import argparse
import json
import sys
from .common import config


class ExactArgumentParser(argparse.ArgumentParser):
    """Subparsers inherit this class: long options must always be exact."""
    def __init__(self, *args, **kwargs):
        kwargs['allow_abbrev'] = False
        super().__init__(*args, **kwargs)


def parser():
    root = ExactArgumentParser(prog='nexvirome2')
    root.add_argument('--version', action='version', version='nexvirome2 0.1.0')
    commands = root.add_subparsers(dest='command', required=True)
    from .research_cli import register
    register(commands)
    commands.add_parser('doctor', help='Check installed external tools and executable startup').add_argument('--extended',action='store_true')
    ref = commands.add_parser('reference').add_subparsers(dest='action', required=True)
    fetch = ref.add_parser('fetch')
    for key in ('accessions', 'config', 'output'):
        fetch.add_argument('--'+key, required=True)
    ref.add_parser('validate').add_argument('--snapshot', required=True)
    panel = commands.add_parser('panel').add_subparsers(dest='action', required=True).add_parser('build')
    panel.add_argument('--reference', required=True)
    panel.add_argument('--metadata', required=True)
    sub = {'panel': panel}
    for name in ('prepare', 'simulate', 'assemble', 'evaluate', 'correct', 'report'):
        sub[name] = commands.add_parser(name)
    for group,actions in {'features':['coverage','propose','neighborhoods'],'latent':['fit'],'orfs':['predict'],
                          'search':['sequence','structure'],'structure':['annotate','database','motif','align'],
                          'bins':['build'],'qc':['run','database']}.items():
        actions_parser = commands.add_parser(group).add_subparsers(dest='action',required=True)
        for action in actions:
            sub[group+'_'+action] = actions_parser.add_parser(action)
    for item in sub.values():
        item.add_argument('--config', required=True)
        item.add_argument('--output', required=True)
        item.add_argument('--threads', type=int)
    sub['prepare'].add_argument('--panel', required=True)
    sub['simulate'].add_argument('--sources', required=True)
    sub['simulate'].add_argument('--seed', type=int)
    sub['assemble'].add_argument('--assembler', choices=['metaspades', 'megahit'], required=True)
    sub['assemble'].add_argument('--memory-gb', type=int)
    for name in ('assemble', 'correct'):
        sub[name].add_argument('--r1', required=True)
        sub[name].add_argument('--r2', required=True)
    for name in ('evaluate', 'correct'):
        sub[name].add_argument('--contigs', required=True)
    sub['evaluate'].add_argument('--sources', required=True)
    sub['evaluate'].add_argument('--alignments', help='Precomputed 11-column BLAST TSV')
    sub['correct'].add_argument('--graph', required=True)
    sub['correct'].add_argument('--paths', required=True)
    sub['correct'].add_argument('--proposals',nargs='*',default=[])
    sub['report'].add_argument('--runs', required=True)
    sub['features_coverage'].add_argument('--contigs',required=True)
    sub['features_coverage'].add_argument('--samples',required=True)
    sub['features_coverage'].add_argument('--regions')
    for key in ('features_propose','latent_fit','bins_build'):
        for field in ('coverage','windows','libraries'):
            sub[key].add_argument('--'+field,required=True)
    for key in ('features_propose','bins_build'):
        sub[key].add_argument('--composition',required=True)
    sub['bins_build'].add_argument('--profiles')
    sub['orfs_predict'].add_argument('--contigs',required=True)
    sub['orfs_predict'].add_argument('--imported',help='Optional trusted ORF coordinate TSV')
    sub['orfs_predict'].add_argument('--regions',help='Optional original-coordinate selection TSV')
    for key in ('search_sequence','search_structure'):
        sub[key].add_argument('--proteins',required=True)
        sub[key].add_argument('--database',required=True)
        sub[key].add_argument('--cache')
    structure_input = sub['search_structure'].add_mutually_exclusive_group(required=True)
    structure_input.add_argument('--model')
    structure_input.add_argument('--structure-map')
    sub['structure_annotate'].add_argument('--orfs',required=True)
    sub['structure_annotate'].add_argument('--metadata',required=True)
    sub['structure_annotate'].add_argument('--sequence-hits')
    sub['structure_annotate'].add_argument('--structure-hits')
    sub['structure_database'].add_argument('--name',choices=['ProstT5','PDB','BFVD'],required=True)
    for key in ('qc_run','qc_database'):
        sub[key].add_argument('--tool',choices=['genomad','checkv'],required=True)
    sub['qc_run'].add_argument('--contigs',required=True)
    sub['qc_run'].add_argument('--database',required=True)
    for field in ('contigs','graph','paths'):
        sub['features_neighborhoods'].add_argument('--'+field,required=True)
    source = sub['structure_motif'].add_mutually_exclusive_group(required=True)
    source.add_argument('--structures')
    source.add_argument('--index')
    sub['structure_motif'].add_argument('--query')
    sub['structure_motif'].add_argument('--residues')
    sub['structure_align'].add_argument('--structures',required=True)
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'doctor':
            from .doctor import check
            result = check(extended=args.extended)
            print(json.dumps(result, indent=2))
            return 0 if result['ok'] else 1
        if args.command == 'reference' and args.action == 'validate':
            from .reference import validate
            print(json.dumps(validate(args.snapshot)))
            return 0
        settings = config(args.config)
        if getattr(args, 'threads', None) is not None:
            if args.threads < 1:
                raise ValueError('threads must be positive')
            settings['threads'] = args.threads
        if getattr(args, 'memory_gb', None) is not None:
            if args.memory_gb < 1:
                raise ValueError('memory-gb must be positive')
            settings['memory_gb'] = args.memory_gb
        from .research_cli import COMMANDS, dispatch
        if args.command in COMMANDS:
            dispatch(args, settings)
        elif args.command == 'reference':
            from .reference import fetch
            fetch(settings, args.accessions, args.output)
        elif args.command == 'panel':
            from .panel import build
            build(settings, args.reference, args.metadata, args.output)
        elif args.command == 'prepare':
            from .benchmark import prepare
            prepare(settings, args.panel, args.output)
        elif args.command == 'simulate':
            from .simulation import simulate
            if args.seed is not None:
                settings['seed'] = args.seed
            simulate(settings, args.sources, args.output)
        elif args.command == 'assemble':
            from .assembly import assemble
            assemble(settings, args.assembler, args.r1, args.r2, args.output)
        elif args.command == 'evaluate':
            from .evaluation import evaluate
            evaluate(settings, args.contigs, args.sources, args.output, args.alignments)
        elif args.command == 'correct':
            from .application import run_correction, CorrectionRequest
            request = CorrectionRequest(args.contigs,args.graph,args.paths,args.r1,args.r2,tuple(args.proposals))
            run_correction(settings,request,args.output)
        elif args.command == 'report':
            from .report import report
            report(settings, args.runs, args.output)
        elif args.command=='features':
            from .features import coverage,propose,neighborhoods
            if args.action=='coverage':
                coverage(settings,args.contigs,args.samples,args.output,args.regions)
            elif args.action=='neighborhoods':
                neighborhoods(settings,args.contigs,args.graph,args.paths,args.output)
            else:
                propose(settings,args.coverage,args.windows,args.libraries,args.composition,args.output)
        elif args.command=='latent':
            from .latent import fit
            fit(settings,args.coverage,args.windows,args.libraries,args.output)
        elif args.command=='orfs':
            from .orfs import predict
            predict(settings,args.contigs,args.output,args.imported,args.regions)
        elif args.command=='search':
            from .application import run_search, SearchRequest
            request = SearchRequest(args.proteins,args.database,'mmseqs' if args.action=='sequence' else 'foldseek',
                                    getattr(args,'model',None),getattr(args,'structure_map',None),args.cache)
            run_search(settings,request,args.output)
        elif args.command=='structure':
            from .searches import download_database,motif,align_structures
            from .structural import annotation
            if args.action=='database':
                download_database(settings,args.name,args.output)
            elif args.action=='motif':
                if bool(args.query)!=bool(args.residues):
                    raise ValueError('--query and --residues must be supplied together')
                motif(settings,args.structures,args.output,args.index,args.query,args.residues)
            elif args.action=='align':
                align_structures(settings,args.structures,args.output)
            else:
                annotation(settings,args.orfs,args.metadata,args.output,args.sequence_hits,args.structure_hits)
        elif args.command=='bins':
            from .binning import build
            build(settings,args.coverage,args.windows,args.libraries,args.composition,args.output,args.profiles)
        elif args.command=='qc':
            from .quality import run,download
            if args.action=='database':
                download(settings,args.tool,args.output)
            else:
                run(settings,args.tool,args.contigs,args.database,args.output)
        return 0
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f'nexvirome2: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
