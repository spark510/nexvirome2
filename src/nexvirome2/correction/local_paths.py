from ..domain.local import LocalPaths

def build_local_paths(graph, candidate, context=600, maximum_paths=128):
    current = candidate['current']
    entrances = graph.incoming(current[1])
    # Include other entrances as competitors so shared upstream sequence
    # cannot falsely assign a fragment from another source to this path.
    all_paths = [current, *candidate['alternatives']]
    for entrance in entrances:
        if entrance != current[0]:
            for exit_node in sorted(graph.links[current[-2]]):
                all_paths.append([entrance, *current[1:-1], exit_node])
    if len(all_paths) > maximum_paths:
        return None
    templates, bounds = {}, {}
    for j, path in enumerate(all_paths):
        seq, offsets = graph.spell(path)
        incoming_overlap = graph.links[path[0]][path[1]]
        left_bound = len(graph.sequence(path[0]))-incoming_overlap
        right_bound = offsets[-1]+graph.links[path[-2]][path[-1]]
        start = max(0, left_bound-context)
        end = min(len(seq), right_bound+context)
        name = f'path_{j}'
        templates[name] = seq[start:end]
        bounds[name] = (left_bound-start, right_bound-start)
    return LocalPaths(templates, bounds, all_paths)
