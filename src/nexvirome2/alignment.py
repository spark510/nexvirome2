"""BLAST with all competing subjects retained and explicit aligned strings."""
from .common import command, version

FIELDS = 'qseqid sseqid pident length qstart qend sstart send bitscore qseq sseq'


def blast(query, subject, output, manifest, name='alignments', threads=1):
    from .common import fasta
    version('blastn', '-version', output, manifest)
    target = output / f'{name}.tsv'
    command(['blastn', '-task', 'blastn', '-query', query, '-subject', subject,
             '-dust', 'no', '-soft_masking', 'false', '-word_size', '11',
             '-max_target_seqs', str(max(1, len(fasta(subject)))),
             '-num_threads', threads, '-outfmt', '6 ' + FIELDS, '-out', target],
            output, manifest, name)
    return target


def read_hits(path):
    hits = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            values = line.rstrip().split('\t')
            if len(values) != 11:
                raise ValueError('Expected BLAST format: ' + FIELDS)
            h = dict(zip(FIELDS.split(), values))
            for key in ('length', 'qstart', 'qend', 'sstart', 'send'):
                h[key] = int(h[key])
            for key in ('pident', 'bitscore'):
                h[key] = float(h[key])
            hits.append(h)
    return hits
