from dataclasses import dataclass


@dataclass(frozen=True)
class CorrectionSettings:
    min_fragments: int = 3
    alternative_fraction: float = .9
    proposal_policy: str = 'annotate'
    proposal_min_score: float = .5
    required_methods: tuple[str, ...] = ()

    def __post_init__(self):
        if self.min_fragments<1 or not 0<self.alternative_fraction<=1:
            raise ValueError('Invalid correction thresholds')
        if self.proposal_policy not in ('annotate','require_any','require_all'):
            raise ValueError('Invalid proposal policy')
        if self.proposal_policy=='require_all' and not self.required_methods:
            raise ValueError('require_all needs explicit required_methods')
        if not 0<=self.proposal_min_score<=1:
            raise ValueError('Invalid proposal minimum score')

    @classmethod
    def from_mapping(cls, settings):
        return cls(settings.get('min_fragments',3),settings.get('alternative_fraction',.9),
                   settings.get('proposal_policy','annotate'),settings.get('proposal_min_score',.5),
                   tuple(settings.get('required_methods',())))
