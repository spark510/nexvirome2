import csv
import gzip
import hashlib
import json
from pathlib import Path

def config(path):
    path = Path(path)
    if path.suffix == '.json':
        return json.loads(path.read_text(encoding='utf-8'))
    import yaml
    return yaml.safe_load(path.read_text(encoding='utf-8')) or {}


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def openseq(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def fasta(path):
    result = {}
    name = None
    with openseq(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
                result[name] = ''
            elif name is None:
                raise ValueError('Sequence before FASTA header')
            else:
                result[name] += line.upper()
    return result


def write_fasta(path, records):
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        for name, seq in records.items():
            handle.write(f'>{name}\n')
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start+80] + '\n')


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN'))[::-1]


def table(path, rows, fields):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def rows(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if fields is None:
            return []
        if len(set(fields)) != len(fields):
            raise ValueError(f'Duplicate TSV column names: {path}')
        if any(not field for field in fields):
            raise ValueError(f'Empty TSV column name: {path}')
        result = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'TSV column count mismatch at line {reader.line_num}: {path}')
            result.append(row)
        return result


def union_length(intervals):
    end = -1
    total = 0
    for left, right in sorted(intervals):
        total += max(0, right - max(left, end))
        end = max(end, right)
    return total


