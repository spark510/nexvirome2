"""Conservative, paired-fragment-supported replacement of bounded graph paths."""
from collections import defaultdict
import math

from .common import fasta, write_fasta, stage, dump, table
from .graph import Graph, contig_paths
from .runtime import CommandRunner
from .correction.support import pair_evidence, validate_score_margin


def bounded_paths(graph, start, end, max_nodes, max_paths, max_states):
    """Enumerate simple paths; any truncation/cycle invalidates this candidate."""
    stack, found, visited = [[start]], [], 0
    while stack:
        path = stack.pop()
        visited += 1
        if visited > max_states:
            return None
        if path[-1] == end:
            found.append(path)
            if len(found) > max_paths:
                return None
            continue
        edges = graph.links.get(path[-1], {})
        if edges and len(path) >= max_nodes:
            return None
        for nxt in sorted(edges, reverse=True):
            if nxt[:-1] in {n[:-1] for n in path}:
                return None
            stack.append(path+[nxt])
    return found


def discover(graph, sequences, paths, max_nodes=8, max_paths=32, max_states=4096):
    candidates, skipped = [], []
    for name, sequence in sequences.items():
        original = paths.get(name)
        try:
            valid = original and None not in original and graph.spell(original)[0] == sequence
        except (KeyError, ValueError):
            valid = False
        if not valid or len({n[:-1] for n in original}) != len(original):
            skipped.append(dict(contig=name, reason='missing_gapped_cyclic_or_mismatched_path'))
            continue
        for left in range(len(original)-1):
            if len(graph.links.get(original[left], {})) < 2:
                continue
            for right in range(left+1, min(len(original), left+max_nodes)):
                alternatives = bounded_paths(graph, original[left], original[right], max_nodes, max_paths, max_states)
                if alternatives is None:
                    skipped.append(dict(contig=name, left=left, right=right, reason='cycle_or_search_limit'))
                    continue
                current = original[left:right+1]
                alternatives = [p for p in alternatives if p != current]
                if not alternatives:
                    continue
                candidates.append(dict(contig=name, left=left, right=right,
                                       original=original, alternatives=alternatives))
                # Nearest fully enumerated reconvergence for this branch.
                break
    return candidates, skipped


def templates_for(graph, sequences, candidate):
    name, original = candidate['contig'], candidate['original']
    left, right = candidate['left'], candidate['right']
    full_paths = [original] + [original[:left]+p+original[right+1:] for p in candidate['alternatives']]
    templates, bounds = {}, {}
    for i, path in enumerate(full_paths):
        sequence, offsets = graph.spell(path)
        exit_index = right if i == 0 else left+len(candidate['alternatives'][i-1])-1
        # Both outside anchors must be spanned, so a fragment connects all
        # changed junctions rather than independently supporting disconnected edges.
        upstream = offsets[left]+len(graph.sequence(path[left]))-graph.links[path[left]][path[left+1]]
        downstream = offsets[exit_index]+graph.links[path[exit_index-1]][path[exit_index]]
        templates[f'path_{i}'] = sequence
        bounds[f'path_{i}'] = (upstream, downstream)
    for i, (other, sequence) in enumerate(sequences.items()):
        if other != name:
            templates[f'background_{i}'] = sequence
            # Background placements compete but cannot support a rewrite.
            bounds[f'background_{i}'] = (0, len(sequence))
    return templates, bounds, full_paths


def coordinate_rows(graph, name, path, original, candidate=None):
    _, offsets = graph.spell(path)
    _, old_offsets = graph.spell(original)
    old_positions = {node:i for i,node in enumerate(original)}
    result = []
    for i, node in enumerate(path):
        overlap = graph.links[path[i-1]][node] if i else 0
        length = len(graph.sequence(node))
        if overlap == length:
            continue
        old_index = i if candidate is None else old_positions.get(node)
        # Retain original coordinates only outside the replaced interior.
        known = old_index is not None and (candidate is None or old_index <= candidate['left'] or old_index >= candidate['right'])
        result.append(dict(output_id=name, output_start=offsets[i]+overlap, output_end=offsets[i]+length,
                           graph_node=node[:-1], strand=node[-1],
                           node_start=overlap if node[-1]=='+' else 0,
                           node_end=length if node[-1]=='+' else length-overlap,
                           original_contig=name if known else None,
                           original_start=old_offsets[old_index]+overlap if known else None,
                           original_end=old_offsets[old_index]+length if known else None,
                           operation='retained' if known else 'graph_replacement'))
    return result


def reconstruct(settings, contigs, graph_file, paths_file, r1, r2, output, runner=None):
    runner = runner or CommandRunner()
    margin = validate_score_margin(settings.get('score_margin', 1))
    minimum, fraction = settings.get('min_fragments', 3), settings.get('alternative_fraction', .9)
    limits = {k:settings.get(k,v) for k,v in dict(max_nodes=8,max_paths=32,max_states=4096,max_fragment_bp=600,threads=4).items()}
    if any(type(v) is not int or v < 1 for v in [minimum,*limits.values()]) or not math.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError('Invalid reconstruction thresholds/limits')
    apply = settings.get('apply', False)
    if type(apply) is not bool:
        raise ValueError('apply must be a boolean')
    with stage(output, settings, [contigs,graph_file,paths_file,r1,r2]) as (out, manifest):
        sequences, graph, paths = fasta(contigs), Graph.load(graph_file), contig_paths(paths_file)
        candidates, skipped = discover(graph,sequences,paths,limits['max_nodes'],limits['max_paths'],limits['max_states'])
        dump(out/'candidates.json', dict(candidates=candidates,skipped=skipped))
        decisions, selected = [], defaultdict(list)
        if candidates:
            runner.version('bowtie2','--version',out,manifest)
        for index, candidate in enumerate(candidates):
            row = dict(candidate=index,contig=candidate['contig'],decision='unresolved')
            decisions.append(row)
            templates, bounds, full_paths = templates_for(graph,sequences,candidate)
            if any(b-a+40 > limits['max_fragment_bp'] for k,(a,b) in bounds.items() if k.startswith('path_')):
                row['reason'] = 'region_exceeds_fragment_span'
                continue
            local = out/f'candidate_{index}'
            local.mkdir()
            write_fasta(local/'templates.fasta',templates)
            dump(local/'paths.json',dict(paths=full_paths,bounds=bounds))
            runner(['bowtie2-build',local/'templates.fasta',local/'index'],out,manifest,f'index_{index}')
            runner(['bowtie2','-x',local/'index','-1',r1,'-2',r2,'-a','--end-to-end',
                    '--very-sensitive','--no-mixed','--no-discordant','-I','0','-X',limits['max_fragment_bp'],
                    '-p',limits['threads'],'-S',local/'reads.sam'],out,manifest,f'mapping_{index}')
            support, ambiguous = pair_evidence(local/'reads.sam',bounds,margin,local/'fragments.json')
            row.update(support=support, ambiguous_fragments=ambiguous)
            ranked = sorted(((support[f'path_{i}'],i) for i in range(1,len(full_paths))),reverse=True)
            count, winner = ranked[0]
            total = sum(support.values())
            unique = len(ranked)==1 or count > ranked[1][0]
            if not unique or count < minimum or not total or count/total < fraction:
                row['reason'] = 'insufficient_unique_alternative_support'
                continue
            if templates[f'path_{winner}'] == sequences[candidate['contig']]:
                row['reason'] = 'sequence_unchanged'
                continue
            if len({n[:-1] for n in full_paths[winner]}) != len(full_paths[winner]):
                row['reason'] = 'replacement_reuses_external_path_nodes'
                continue
            row.update(decision='proposed',reason='spanning_fragment_support',selected_path=full_paths[winner])
            selected[candidate['contig']].append((candidate,full_paths[winner],row))
        result, mapping, proposed_sequences, output_paths = dict(sequences), [], {}, dict(paths)
        for name, choices in selected.items():
            if len(choices) != 1:
                for _,_,row in choices:
                    row.update(decision='unresolved',reason='multiple_rewrites_require_joint_validation')
                continue
            candidate, path, row = choices[0]
            proposed_sequences[name] = graph.spell(path)[0]
            if apply:
                result[name], output_paths[name] = proposed_sequences[name], path
                row['decision'] = 'reconstruct'
        if any(r['decision']=='reconstruct' for r in decisions):
            # Recheck the proposed assembly as a whole: competing rewritten
            # contigs can invalidate support that was unique in a local run.
            verification = out/'verification'
            verification.mkdir()
            write_fasta(verification/'provisional.fasta',result)
            bounds = {name:(0,len(seq)) for name,seq in result.items()}
            for name, choices in selected.items():
                if len(choices)==1 and choices[0][2]['decision']=='reconstruct':
                    candidate,path,_ = choices[0]
                    _, candidate_bounds, full_paths = templates_for(graph,sequences,candidate)
                    bounds[name] = candidate_bounds[f'path_{full_paths.index(path)}']
            runner(['bowtie2-build',verification/'provisional.fasta',verification/'index'],out,manifest,'verification_index')
            runner(['bowtie2','-x',verification/'index','-1',r1,'-2',r2,'-a','--end-to-end',
                    '--very-sensitive','--no-mixed','--no-discordant','-I','0','-X',limits['max_fragment_bp'],
                    '-p',limits['threads'],'-S',verification/'reads.sam'],out,manifest,'verification_mapping')
            verified, ambiguous = pair_evidence(verification/'reads.sam',bounds,margin,verification/'fragments.json')
            dump(verification/'support.json',dict(support=verified,ambiguous=ambiguous))
            failed = [name for name, choices in selected.items()
                      if len(choices)==1 and choices[0][2]['decision']=='reconstruct' and verified[name]<minimum]
            # Rolling back one contig could introduce a competitor for another.
            # Withhold the whole provisional batch instead of assuming safety.
            for name, choices in selected.items():
                if len(choices)==1 and choices[0][2]['decision']=='reconstruct':
                    row=choices[0][2]
                    row['verification_support']=verified[name]
                    if failed:
                        result[name],output_paths[name]=sequences[name],paths[name]
                        row.update(decision='unresolved',reason='batch_remapping_failed',failed_contigs=failed)
                    else:
                        row['reason']='spanning_support_and_batch_remapping'
        path_validation = {}
        for name in sequences:
            choices = selected.get(name,[])
            applied = len(choices)==1 and choices[0][2]['decision']=='reconstruct'
            path = output_paths.get(name)
            if path and None not in path:
                try:
                    if graph.spell(path)[0] == result[name]:
                        mapping.extend(coordinate_rows(graph,name,path,paths[name],choices[0][0] if applied else None))
                        path_validation[name] = dict(valid=True,reason='sequence_matches_output')
                        continue
                except (KeyError, ValueError):
                    pass
            if applied:
                raise ValueError(f'Reconstructed path does not match output sequence: {name}')
            output_paths[name] = None
            path_validation[name] = dict(valid=False,reason='missing_gapped_invalid_or_mismatched_path')
            mapping.append(dict(output_id=name,output_start=0,output_end=len(result[name]),
                                original_contig=name,original_start=0,original_end=len(result[name]),operation='retained_unmapped'))
        write_fasta(out/'reconstructed.fasta',result)
        write_fasta(out/'proposed.fasta',proposed_sequences)
        dump(out/'output_paths.json',{n:output_paths.get(n) for n in sequences})
        dump(out/'output_path_validation.json',path_validation)
        dump(out/'decisions.json',decisions)
        table(out/'graph_coordinate_map.tsv',mapping,['output_id','output_start','output_end','graph_node','node_start','node_end','strand','original_contig','original_start','original_end','operation'])
        dump(out/'summary.json',dict(applied=sum(r['decision']=='reconstruct' for r in decisions),
             proposed=len(proposed_sequences),truth_used=False,mode='apply' if apply else 'report_only',
             limitation='One bounded rewrite per contig; PE evidence only; no long-range phasing or lost-path recovery'))
