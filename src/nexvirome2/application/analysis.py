from dataclasses import dataclass


@dataclass(frozen=True)
class CorrectionRequest:
    """Inference inputs only: simulation sources and truth are not accepted."""
    contigs: str
    graph: str
    paths: str
    r1: str
    r2: str
    proposals: tuple[str, ...] = ()


@dataclass(frozen=True)
class SearchRequest:
    proteins: str
    database: str
    backend: str
    model: str | None = None
    structure_map: str | None = None
    cache: str | None = None


def run_correction(settings, request, output, service=None):
    from ..correction.service import CorrectionService
    service = service or CorrectionService()
    return service.run(settings,request.contigs,request.graph,request.paths,request.r1,request.r2,output,request.proposals)


def run_search(settings, request, output, service=None):
    from ..search.service import SearchService
    from ..search.backends import backend_for
    service = service or SearchService(backend_for(request.backend))
    if service.backend.name != request.backend:
        raise ValueError('Injected search backend does not match request')
    return service.run(settings,request.proteins,request.database,output,request.model,request.structure_map,request.cache)


def ablation_settings(settings, variant, threads):
    if variant not in ('baseline','simple','nmf','structure','combined','multiview'):
        raise ValueError('Unknown correction ablation')
    result = dict(settings,threads=threads,proposal_policy='annotate' if variant=='baseline' else 'require_any')
    if variant=='combined':
        result.update(proposal_policy='require_all',required_methods=['nmf','structural_lineage'])
    return result
