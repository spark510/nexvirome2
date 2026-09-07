"""Accession-exact reference acquisition; no biological filtering at download."""
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from ..common import dump, fasta, rows, sha, stage, table, write_fasta


def accessions(path):
    result = []
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        accession = line.strip()
        if not accession:
            continue
        if not re.fullmatch(r'(?:NC|AC|NZ|NM|NR|XM|XR)_\d+\.\d+', accession):
            raise ValueError(f'Expected one versioned RefSeq nucleotide accession per line: {accession}')
        if accession not in result:
            result.append(accession)
    if not result:
        raise ValueError('Empty accession list')
    return result


def download(url, target, payload=None):
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, data=payload, headers={'User-Agent': 'NexVirome2/0.1'})
            with urllib.request.urlopen(request, timeout=180) as response:
                content = response.read()
            Path(target).write_bytes(content)
            return
        except (OSError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def parse_genbank(path):
    records, metadata = {}, []
    root = ET.parse(path).getroot()
    for record in root.findall('.//GBSeq'):
        accession = record.findtext('GBSeq_accession-version')
        seq = (record.findtext('GBSeq_sequence') or '').upper()
        if not accession or not seq:
            raise ValueError('GenBank record lacks accession or sequence')
        if accession in records:
            raise ValueError(f'Duplicate accession: {accession}')
        records[accession] = seq
        qualifiers = {}
        for feature in record.findall('.//GBFeature'):
            if feature.findtext('GBFeature_key') == 'source':
                for qualifier in feature.findall('.//GBQualifier'):
                    qualifiers.setdefault(qualifier.findtext('GBQualifier_name'), []).append(qualifier.findtext('GBQualifier_value') or '')
        definition = record.findtext('GBSeq_definition') or ''
        taxa = [v[6:] for v in qualifiers.get('db_xref', []) if v.startswith('taxon:')]
        completeness = 'complete' if re.search(r'complete (?:genome|sequence)', definition, re.I) else 'unknown'
        metadata.append({'accession': accession, 'length': len(seq),
                         'organism': record.findtext('GBSeq_organism'), 'taxonomy': record.findtext('GBSeq_taxonomy'),
                         'taxid': ';'.join(taxa), 'molecule': record.findtext('GBSeq_moltype'),
                         'topology': record.findtext('GBSeq_topology'), 'segment': ';'.join(qualifiers.get('segment', [])),
                         'completeness': completeness, 'definition': definition})
    return records, metadata


def fetch(settings, accession_file, output):
    requested = accessions(accession_file)
    with stage(output, settings, [accession_file]) as (out, manifest):
        raw = out / 'raw'
        raw.mkdir()
        (out / 'accessions.txt').write_text('\n'.join(requested) + '\n')
        sequences, metadata = {}, []
        endpoint = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        for batch, start in enumerate(range(0, len(requested), 100)):
            fields = {'db': 'nuccore', 'id': ','.join(requested[start:start+100]),
                      'rettype': 'gb', 'retmode': 'xml', 'tool': 'nexvirome2'}
            if settings.get('email'):
                fields['email'] = settings['email']
            path = raw / f'batch_{batch:05}.xml'
            download(endpoint, path, urllib.parse.urlencode(fields).encode())
            batch_seq, batch_meta = parse_genbank(path)
            if sequences.keys() & batch_seq.keys():
                raise ValueError('Duplicate accession across API batches')
            sequences.update(batch_seq)
            metadata.extend(batch_meta)
            time.sleep(0.4)
        reconciliation = {'missing': sorted(set(requested)-sequences.keys()),
                          'extra': sorted(sequences.keys()-set(requested))}
        dump(out / 'reconciliation.json', reconciliation)
        write_fasta(out / 'reference.fasta', sequences)
        fields = ['accession', 'length', 'organism', 'taxonomy', 'taxid', 'molecule', 'topology', 'segment', 'completeness', 'definition']
        table(out / 'metadata.tsv', metadata, fields)
        manifest.update(source_query=settings.get('source_query'), acquisition_endpoint=endpoint,
                        accession_count=len(sequences), completeness_method='GenBank definition; unknown is preserved')
        if reconciliation['missing'] or reconciliation['extra']:
            raise ValueError('Accession/version mismatch; see reconciliation.json')


def validate(snapshot):
    import json
    snapshot = Path(snapshot)
    manifest = json.loads((snapshot / 'manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Snapshot is not complete')
    for name, checksum in manifest['outputs'].items():
        if sha(snapshot / name) != checksum:
            raise ValueError(f'Checksum mismatch: {name}')
    expected = set(accessions(snapshot / 'accessions.txt'))
    sequences = fasta(snapshot / 'reference.fasta')
    metadata = rows(snapshot / 'metadata.tsv')
    if expected != set(sequences) or expected != {r['accession'] for r in metadata} or len(metadata) != len(expected):
        raise ValueError('Accession sets differ')
    for row in metadata:
        if len(sequences[row['accession']]) != int(row['length']):
            raise ValueError('Metadata length mismatch')
    return {'valid': True, 'accessions': len(expected)}
