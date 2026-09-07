"""Exact shared-kmer candidates, disk-backed to bound index memory."""
import itertools
import math
import sqlite3
from collections import Counter
from pathlib import Path

from .alignment import blast, read_hits
from .common import dump, fasta, revcomp, rows, stage, table, write_fasta


def low_complexity(seq):
    counts = Counter(seq)
    entropy = -sum((n/len(seq))*math.log2(n/len(seq)) for n in counts.values())
    return entropy < 1.5


def kmers(seq, k):
    counts = Counter()
    for i in range(len(seq)-k+1):
        word = seq[i:i+k]
        if set(word) <= set('ACGT'):
            counts[min(word, revcomp(word))] += 1
    return counts


def build(settings, reference, metadata, output):
    with stage(output, settings, [reference, metadata]) as (out, manifest):
        all_seq = fasta(reference)
        meta = {r['accession']: r for r in rows(metadata)}
        selected = {a: seq for a, seq in all_seq.items() if len(seq) >= settings.get('min_length_bp', 1000)
                    and a not in set(settings.get('exclude_accessions', []))
                    and (not settings.get('complete_only', False) or meta.get(a, {}).get('completeness') == 'complete')}
        k = settings.get('k', 31)
        db = sqlite3.connect(out / 'kmers.sqlite')
        try:
            db.execute('CREATE TABLE kmers (word TEXT, accession TEXT, repeated INTEGER, low INTEGER)')
            for accession, seq in sorted(selected.items()):
                db.executemany('INSERT INTO kmers VALUES (?,?,?,?)',
                               ((word, accession, n > 1, low_complexity(word)) for word, n in kmers(seq, k).items()))
                db.commit()
            db.execute('CREATE INDEX word_idx ON kmers(word)')
            query = '''SELECT a.accession,b.accession,count(*),
                       sum(CASE WHEN a.low=1 OR b.low=1 THEN 1 ELSE 0 END),
                       sum(CASE WHEN a.repeated=1 OR b.repeated=1 THEN 1 ELSE 0 END)
                       FROM kmers a JOIN kmers b ON a.word=b.word AND a.accession<b.accession
                       GROUP BY a.accession,b.accession ORDER BY count(*) DESC,a.accession,b.accession'''
            candidates = [{'a': a, 'b': b, 'shared_kmers': n, 'low_complexity_kmers': low, 'repeated_kmers': rep}
                          for a, b, n, low, rep in db.execute(query)]
        finally:
            db.close()
        table(out / 'candidate_pairs.tsv', candidates, ['a', 'b', 'shared_kmers', 'low_complexity_kmers', 'repeated_kmers'])
        candidates.sort(key=lambda c: (-(c['shared_kmers']-c['low_complexity_kmers']), c['a'], c['b']))
        accepted = []
        for index, candidate in enumerate(candidates):
            if candidate['shared_kmers'] == candidate['low_complexity_kmers']:
                continue
            a, b = candidate['a'], candidate['b']
            query_fa, subject_fa = out / 'query.fasta', out / 'subject.fasta'
            write_fasta(query_fa, {a: selected[a]})
            write_fasta(subject_fa, {b: selected[b]})
            hits = read_hits(blast(query_fa, subject_fa, out, manifest, f'pair_{index}', settings.get('threads', 1)))
            flank = settings.get('flank_bp', 150)
            # Each side must contain a kmer absent anywhere in the competing genome.
            ka, kb = kmers(selected[a], k), kmers(selected[b], k)
            useful = []
            for hit in hits:
                ql, qr = sorted((hit['qstart']-1, hit['qend']-1))
                sl, sr = sorted((hit['sstart']-1, hit['send']-1))
                if min(ql, sl, len(selected[a])-qr-1, len(selected[b])-sr-1) < flank:
                    continue
                flanks = [(selected[a][ql-flank:ql], kb), (selected[a][qr+1:qr+1+flank], kb),
                          (selected[b][sl-flank:sl], ka), (selected[b][sr+1:sr+1+flank], ka)]
                if all(set(kmers(seq, k))-set(other) for seq, other in flanks):
                    useful.append(hit)
            if useful:
                accepted.append({**candidate, 'shared_regions': useful})
            if len(accepted) == settings.get('high_risk_pairs', 2):
                break
        if len(accepted) < settings.get('high_risk_pairs', 2):
            raise ValueError('Too few naturally shared regions with discriminating flanks; inspect candidate_pairs.tsv')
        used = {x for pair in accepted for x in (pair['a'], pair['b'])}
        shared = {(r['a'], r['b']): r['shared_kmers'] for r in candidates}
        controls = itertools.combinations(sorted(set(selected)-used), 2)
        control = min(controls, key=lambda p: (shared.get(p, 0), p), default=None)
        if control is None:
            raise ValueError('Need two additional references for an independent low-sharing control')
        accepted.append({'a': control[0], 'b': control[1], 'shared_kmers': shared.get(control, 0), 'control': True})
        panel_dir = out / 'sources'
        panel_dir.mkdir()
        for accession in sorted({x for p in accepted for x in (p['a'], p['b'])}):
            write_fasta(panel_dir / f'{accession}.fasta', {accession: selected[accession]})
        dump(out / 'panel.json', {'pairs': accepted, 'metadata': {a: meta.get(a, {}) for a in selected if a in {x for p in accepted for x in (p['a'],p['b'])}}})
