import random
from collections import defaultdict
from ..common import dump, fasta, stage, table, write_fasta
from ..graph import Graph, candidates, contig_paths
from ..evidence.proposals import proposals
from ..runtime import CommandRunner
from .support import pair_evidence, validate_score_margin
from .splitting import split_sequences
from .policy import EvidencePolicy, CorrectionPolicy
from .local_paths import build_local_paths
from ..domain import PairSupport

class CorrectionService:
    """Compose graph mechanics, an execution dependency, and a decision policy."""
    def __init__(self, runner=None, version=None, policy: CorrectionPolicy | None = None):
        runner = runner or CommandRunner()
        self.command = runner
        self.version = version or runner.version
        self.policy = policy

    def run(self, settings, contigs, graph_file, paths_file, r1, r2, output, proposal_files=()):
        command, version = self.command, self.version
        policy = self.policy or EvidencePolicy(settings)
        validate_score_margin(settings.get('score_margin', 1))
        with stage(output, settings, [contigs, graph_file, paths_file, r1, r2, *proposal_files]) as (out, manifest):
            sequences = fasta(contigs)
            extra = proposals(proposal_files,sequences)
            graph = Graph.load(graph_file)
            proposed = candidates(graph, sequences, contig_paths(paths_file), settings.get('max_repeat_nodes', 20))
            version('bowtie2', '--version', out, manifest)
            version('samtools', '--version', out, manifest)
            cuts, evidence = defaultdict(list), []
            context = settings.get('anchor_context_bp', 600)
            maximum = settings.get('max_fragment_bp', 600)
            for i, candidate in enumerate(proposed):
                row = {'candidate': i, **candidate, 'current_support': 0, 'alternative_support': 0,
                       'ambiguous_pairs': 0, 'decision': 'unresolved'}
                if not candidate['usable']:
                    evidence.append(row)
                    continue
                current = candidate['current']
                local_paths = build_local_paths(graph,candidate,context,settings.get('max_local_paths',128))
                if local_paths is None:
                    row['reason'] = 'too_many_competing_paths'
                    evidence.append(row)
                    continue
                templates, bounds, all_paths = local_paths.templates, local_paths.bounds, local_paths.paths
                if bounds['path_0'][1]-bounds['path_0'][0]+40 > maximum:
                    row['reason'] = 'repeat_exceeds_fragment_span'
                    evidence.append(row)
                    continue
                local = out / f'candidate_{i}'
                local.mkdir()
                targets = local / 'templates.fasta'
                write_fasta(targets, templates)
                dump(local / 'paths.json', {'paths': all_paths, 'bounds': bounds})
                index = local / 'index'
                command(['bowtie2-build', targets, index], out, manifest, f'index_{i}')
                sam = local / 'alignments.sam'
                command(['bowtie2', '-x', index, '-1', r1, '-2', r2, '-a', '--end-to-end',
                         '--very-sensitive', '--no-mixed', '--no-discordant', '-I', '0', '-X', maximum,
                         '-p', settings.get('threads', 4), '-S', sam], out, manifest, f'mapping_{i}')
                support, ambiguous = pair_evidence(sam, bounds, settings.get('score_margin', 1), local/'fragment_candidates.json')
                dump(local/'path_support.json', [{'path':name,'independent_fragments':count,
                     'nodes':all_paths[int(name.split('_')[-1])]} for name,count in
                     sorted(support.items(),key=lambda item:(-item[1],item[0]))])
                current_support = support['path_0']
                # Require a single alternative to meet the independent fragment threshold.
                alternative_support = max((support[f'path_{j}'] for j in range(1, 1+len(candidate['alternatives']))), default=0)
                outcome = policy.decide(PairSupport(current_support, alternative_support), candidate, extra)
                read_call, call, supporting_methods = outcome.read_decision, outcome.decision, list(outcome.methods)
                row.update(read_decision=read_call,proposal_methods=supporting_methods)
                command(['samtools', 'sort', '-o', local / 'alignments.bam', sam], out, manifest, f'sort_{i}')
                command(['samtools', 'depth', '-a', local / 'alignments.bam'], out, manifest, f'depth_{i}', stdout=local / 'depth.tsv')
                row.update(current_support=current_support, alternative_support=alternative_support,
                           ambiguous_pairs=ambiguous, decision=call,
                           reason='alternative_fragment_support' if call == 'split' else 'current_support' if call == 'retain' else
                           'proposal_gate_not_met' if read_call=='split' else 'insufficient_discriminating_evidence')
                row['coverage'] = {n: graph.coverage.get(n[:-1]) for n in current}
                if call == 'split':
                    cuts[candidate['contig']].append(candidate['cut'])
                evidence.append(row)
            corrected, mapping = split_sequences(sequences, cuts)
            write_fasta(out / 'corrected.fasta', corrected)
            table(out / 'coordinate_map.tsv', mapping, ['fragment', 'original', 'start', 'end'])
            table(out / 'evidence.tsv', evidence, ['candidate', 'contig', 'cut', 'repeat_length', 'current_support', 'alternative_support', 'ambiguous_pairs', 'read_decision', 'proposal_methods', 'decision', 'reason'])
            dump(out / 'evidence.json', evidence)
            matched = [p for p in extra if any(r.get('usable') and r['contig']==p['contig'] and p['start']<=r.get('cut',-1)<=p['end'] for r in evidence)]
            dump(out/'proposal_audit.json',{'proposals':extra,'matched':len(matched),'unmatched':len(extra)-len(matched),
                                           'policy':settings.get('proposal_policy','annotate'),
                                           'unmatched_action':'unresolved; no forced split without a validated graph path'})
            rng = random.Random(settings.get('random_control_seed', 20260905))
            random_cuts = {name: rng.sample(range(1, len(sequences[name])), len(set(points))) for name, points in cuts.items()}
            control, control_map = split_sequences(sequences, random_cuts)
            write_fasta(out / 'random_control.fasta', control)
            table(out / 'random_coordinate_map.tsv', control_map, ['fragment', 'original', 'start', 'end'])
            dump(out / 'summary.json', {'split_count': sum(len(set(p)) for p in cuts.values()),
                                        'unresolved': sum(r['decision']=='unresolved' for r in evidence),
                                        'truth_used': False})


