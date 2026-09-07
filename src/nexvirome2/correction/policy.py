from dataclasses import dataclass
from typing import Protocol
from ..domain import BreakpointProposal, PairSupport, CorrectionDecision
from ..domain.settings import CorrectionSettings

def decision(current, alternative, minimum=3, fraction=0.9):
    if minimum < 1 or not 0 < fraction <= 1 or min(current, alternative) < 0:
        raise ValueError('Invalid correction thresholds')
    total = current+alternative
    if alternative >= minimum and total and alternative/total >= fraction:
        return 'split'
    if current >= minimum:
        return 'retain'
    return 'unresolved'


def gate(read_decision,candidate,extra,settings):
    options = settings if isinstance(settings,CorrectionSettings) else CorrectionSettings.from_mapping(settings)
    proposals = [p if isinstance(p,BreakpointProposal) else BreakpointProposal.from_record(p) for p in extra]
    methods = sorted({p.method for p in proposals if p.matches(candidate['contig'],candidate.get('cut',-1),options.proposal_min_score)})
    policy = options.proposal_policy
    required = set(options.required_methods)
    supported = bool(methods) if policy=='require_any' else required <= set(methods)
    if read_decision=='split' and policy!='annotate' and not supported:
        return 'unresolved',methods
    return read_decision,methods


class CorrectionPolicy(Protocol):
    def decide(self, support: PairSupport, candidate, proposals) -> CorrectionDecision: ...


@dataclass(frozen=True)
class EvidencePolicy:
    settings: CorrectionSettings

    def __post_init__(self):
        if not isinstance(self.settings,CorrectionSettings):
            object.__setattr__(self,'settings',CorrectionSettings.from_mapping(self.settings))

    def decide(self, support, candidate, proposals):
        read_call = decision(support.current, support.alternative, self.settings.min_fragments,
                             self.settings.alternative_fraction)
        call, methods = gate(read_call, candidate, proposals, self.settings)
        return CorrectionDecision(read_call, call, tuple(methods))
