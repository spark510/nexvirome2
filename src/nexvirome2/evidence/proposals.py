from ..io.formats import rows
from ..domain import BreakpointProposal

def proposals(paths,sequences):
    result = []
    for path in paths:
        for row in rows(path):
            proposal = BreakpointProposal.from_record(row)
            if proposal.contig not in sequences or proposal.end > len(sequences[proposal.contig]):
                raise ValueError('Invalid candidate interval or confidence')
            result.append({**row, **proposal.to_record()})
    return result


