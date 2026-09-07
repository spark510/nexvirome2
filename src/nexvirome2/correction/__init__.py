"""Compatibility entry points for correction."""
from ..common import command, version
from .support import cigar_length, pair_evidence
from .splitting import split_sequences
from .policy import decision
from .service import CorrectionService

def correct(settings, contigs, graph_file, paths_file, r1, r2, output, proposal_files=()):
    return CorrectionService(command, version).run(settings, contigs, graph_file, paths_file, r1, r2, output, proposal_files)
