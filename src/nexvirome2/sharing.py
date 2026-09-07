"""Position-level exact nucleotide sharing atlas backed by an on-disk seed index."""
from collections import Counter
from itertools import groupby
import math
import sqlite3
from .common import fasta, revcomp, stage, table, sha, dump, rows
from .taxonomy import TaxonomyResolver


def seeds(sequence, k, circular=False):
    if circular and len(sequence)>=k:
        sequence += sequence[:k-1]
    for position in range(max(0,len(sequence)-k+1)):
        word = sequence[position:position+k]
        if set(word) <= set('ACGT'):
            yield position, min(word,revcomp(word))


def low_complexity(word, threshold=1.5):
    return -sum((n/len(word))*math.log2(n/len(word)) for n in Counter(word).values()) < threshold


def build(settings, reference, taxonomy, output, metadata=None):
    k, rank = int(settings.get('k',31)), settings.get('rank','species')
    if k < 3:
        raise ValueError('k must be at least 3')
    with stage(output, settings, [reference,taxonomy]+([metadata] if metadata else [])) as (out,_):
        sequences = fasta(reference)
        topology={}
        if metadata:
            for row in rows(metadata):
                if row['accession'] in topology or row['accession'] not in sequences:
                    raise ValueError('Duplicate/unknown topology accession')
                if row['topology'] not in ('linear','circular'):
                    raise ValueError('Topology must be explicitly linear or circular')
                topology[row['accession']]=row['topology']
        resolver = TaxonomyResolver.load(taxonomy)
        taxa = {acc:resolver.taxon(acc,rank) for acc in sequences}
        digest = sha(reference)
        connection = sqlite3.connect(out/'seeds.sqlite')
        try:
            connection.execute('CREATE TABLE occurrence (seed TEXT, accession TEXT, position INTEGER)')
            for acc,seq in sequences.items():
                connection.executemany('INSERT INTO occurrence VALUES (?,?,?)',
                                       ((seed,acc,pos) for pos,seed in seeds(seq,k,topology.get(acc)=='circular')))
            connection.execute('CREATE INDEX seed_index ON occurrence(seed)')
            connection.commit()
            counts = Counter()
            def evidence():
                cursor = connection.execute('SELECT seed,accession,position FROM occurrence ORDER BY seed,accession,position')
                for seed,group in groupby(cursor, key=lambda r:r[0]):
                    occurrences = list(group)
                    accessions = sorted({r[1] for r in occurrences})
                    distinct = sorted({taxa[a] for a in accessions if taxa[a] is not None})
                    unknown = any(taxa[a] is None for a in accessions)
                    multiplicity = Counter(r[1] for r in occurrences)
                    for _,acc,pos in occurrences:
                        counts['seed_positions'] += 1
                        counts['cross_taxon_positions'] += int(len(distinct)>1)
                        record = dict(accession=acc,start=pos,end=min(pos+k,len(sequences[acc])),rank=rank,taxa_count=len(distinct),
                                   taxa=';'.join(distinct),accessions=';'.join(accessions),
                                   unknown_taxonomy=str(unknown).lower(),
                                   repeat=str(multiplicity[acc]>1).lower(),
                                   low_complexity=str(low_complexity(seed)).lower(),
                                   evidence_id=f'{acc}:{pos}:{k}',reference_sha256=digest,method='exact_canonical_kmer',k=k)
                        yield record
                        if pos+k>len(sequences[acc]):
                            yield {**record,'start':0,'end':pos+k-len(sequences[acc])}
            table(out/'atlas.tsv', evidence(), ['accession','start','end','rank','taxa_count','taxa',
                  'accessions','unknown_taxonomy','repeat','low_complexity','evidence_id','reference_sha256','method','k'])
            dump(out/'summary.json',dict(counts,method='exact_canonical_kmer',k=k,rank=rank,
                 taxonomy_sha256=sha(taxonomy),reference_sha256=digest,
                 circular_accessions=sorted(a for a in topology if topology[a]=='circular'),
                 limitations=['Exact seeds miss divergent homologs','Circularity requires explicit metadata; records shorter than k have no seeds',
                              'Seed multiplicity is not empirical read discriminability']))
        finally:
            connection.close()
