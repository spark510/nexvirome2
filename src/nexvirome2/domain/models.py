from dataclasses import dataclass, asdict
import math


@dataclass(frozen=True)
class SequenceInterval:
    """Zero-based half-open sequence interval [start, end)."""
    start: int
    end: int

    def __post_init__(self):
        if not 0 <= self.start < self.end:
            raise ValueError('Invalid sequence interval')

    def overlaps(self, other):
        return self.start < other.end and self.end > other.start


@dataclass(frozen=True)
class CutInterval:
    """Inclusive possible cut positions; start == end denotes an exact cut."""
    start: int
    end: int

    def __post_init__(self):
        if not 0 <= self.start <= self.end:
            raise ValueError('Invalid candidate interval')

    def contains(self, cut):
        return self.start <= cut <= self.end


@dataclass(frozen=True)
class BreakpointProposal:
    contig: str
    start: int
    end: int
    method: str
    score: float

    def __post_init__(self):
        CutInterval(self.start, self.end)
        if not self.contig or not self.method:
            raise ValueError('Candidate contig and method are required')
        if not math.isfinite(self.score) or not 0 <= self.score <= 1:
            raise ValueError('Invalid candidate confidence')

    @classmethod
    def from_record(cls, row):
        return cls(row['contig'], int(row['start']), int(row['end']), row['method'], float(row['score']))

    def matches(self, contig, cut, minimum):
        return self.contig == contig and CutInterval(self.start, self.end).contains(cut) and self.score >= minimum

    def to_record(self):
        return asdict(self)


@dataclass(frozen=True)
class SearchHit:
    """Search identity and coverage fractions; separate from BLAST percentage identity."""
    query: str
    target: str
    fident: float
    qcov: float
    tcov: float
    evalue: float
    bits: float

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.fident,self.qcov,self.tcov,self.evalue,self.bits)) or min(self.evalue,self.bits)<0:
            raise ValueError('Invalid search score')
        if any(not 0 <= v <= 1 for v in (self.fident,self.qcov,self.tcov)):
            raise ValueError('fident/qcov/tcov must be fractions, not percentages')

    def to_record(self):
        return asdict(self)


@dataclass(frozen=True)
class PairSupport:
    current: int
    alternative: int

    def __post_init__(self):
        if min(self.current,self.alternative)<0:
            raise ValueError('Invalid correction thresholds')


@dataclass(frozen=True)
class CorrectionDecision:
    read_decision: str
    decision: str
    methods: tuple[str, ...] = ()

    def __post_init__(self):
        if self.read_decision not in {'split','retain','unresolved'} or self.decision not in {'split','retain','unresolved'}:
            raise ValueError('Invalid correction decision')
        if self.decision=='split' and self.read_decision!='split':
            raise ValueError('Proposal evidence cannot create a read-unsupported split')
