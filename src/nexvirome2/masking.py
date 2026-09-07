"""Reference-only, coordinate-preserving masks. Original FASTA is immutable."""
from dataclasses import dataclass
import json
import random
from .common import fasta, openseq, write_fasta, rows, table, sha, stage, dump


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


@dataclass(frozen=True)
class MaskSet:
    reference_sha256: str
    lengths: dict
    intervals: tuple
    policy: dict

    @classmethod
    def create(cls, reference, records, policy):
        sequences = fasta(reference)
        if not sequences or any(not seq for seq in sequences.values()):
            raise ValueError('Empty reference sequence')
        if not policy.get('version') or not policy.get('rank') or not isinstance(policy.get('conditions'), dict):
            raise ValueError('Mask policy requires version, rank and conditions')
        validated = []
        for index, row in enumerate(records):
            accession = row['accession']
            start, end = int(row['start']), int(row['end'])
            if accession not in sequences or not 0 <= start < end <= len(sequences[accession]):
                raise ValueError(f'Invalid original-coordinate mask: {accession}:{start}-{end}')
            validated.append(dict(accession=accession, start=start, end=end,
                                  event_id=str(row.get('event_id') or f'm{index}'),
                                  evidence=str(row.get('evidence') or 'imported')))
        return cls(sha(reference), {k: len(v) for k,v in sequences.items()}, tuple(validated), dict(policy))

    @classmethod
    def load(cls, path, reference):
        payload = json.loads(open(path, encoding='utf-8').read())
        if payload.get('schema') != 1 or payload.get('coordinate_system') != '0-based-half-open':
            raise ValueError('Unsupported mask schema/coordinates')
        result = cls.create(reference, payload['intervals'], payload['policy'])
        if result.reference_sha256 != payload['reference_sha256'] or result.lengths != payload['lengths']:
            raise ValueError('Mask reference checksum/length mismatch')
        return result

    def save(self, path):
        dump(path, dict(schema=1, coordinate_system='0-based-half-open',
                        reference_sha256=self.reference_sha256, lengths=self.lengths,
                        intervals=self.intervals, policy=self.policy))

    def merged(self, accession):
        return merge_intervals((r['start'],r['end']) for r in self.intervals if r['accession']==accession)

    def render(self, reference, mode):
        if sha(reference) != self.reference_sha256:
            raise ValueError('Mask reference checksum mismatch')
        if mode not in ('hard', 'soft', 'unmasked'):
            raise ValueError('Mask mode must be hard, soft or unmasked')
        sequences = fasta(reference)
        if mode != 'unmasked':
            for accession, seq in sequences.items():
                chars = list(seq)
                for start, end in self.merged(accession):
                    chars[start:end] = 'N'*(end-start) if mode=='hard' else seq[start:end].lower()
                sequences[accession] = ''.join(chars)
        return sequences


def _case_fasta(path):
    result, name = {}, None
    with openseq(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError('Duplicate masked FASTA ID')
                result[name] = ''
            elif name is None:
                raise ValueError('Sequence before header')
            else:
                result[name] += line
    return result


def import_mask(settings, reference, source, output):
    with stage(output, settings, [reference, source]) as (out, _):
        fmt = settings.get('format', 'intervals')
        if fmt == 'intervals':
            records = rows(source)
        elif fmt in ('hard', 'soft'):
            original, masked = fasta(reference), _case_fasta(source)
            if original.keys() != masked.keys():
                raise ValueError('Masked FASTA accessions differ from original')
            records = []
            for accession, seq in original.items():
                view = masked[accession]
                if len(view) != len(seq):
                    raise ValueError('Masking must preserve sequence length')
                positions = []
                for i, (a,b) in enumerate(zip(seq,view)):
                    if fmt=='hard':
                        if b.upper() != a and b != 'N':
                            raise ValueError('Masked FASTA alters an unmasked nucleotide')
                        is_masked = b=='N' and a!='N'
                    else:
                        if b.upper()!=a:
                            raise ValueError('Soft mask alters nucleotide sequence')
                        is_masked = b.islower()
                    if is_masked:
                        positions.append((i,i+1))
                records.extend(dict(accession=accession,start=a,end=b,evidence=f'{fmt}-fasta')
                               for a,b in merge_intervals(positions))
        else:
            raise ValueError('Supported mask imports: intervals, hard, soft')
        mask = MaskSet.create(reference, records, settings['policy'])
        mask.save(out/'mask.json')
        dump(out/'audit.json', {'format':fmt, 'original_preserved':True,
             'native_N_mask_status':'unknown' if fmt=='hard' else 'not_applicable',
             'masked_bases':sum(b-a for acc in mask.lengths for a,b in mask.merged(acc))})


def export_mask(settings, reference, mask_file, output):
    with stage(output, settings, [reference, mask_file]) as (out, _):
        mask = MaskSet.load(mask_file, reference)
        mask.save(out/'mask.json')
        table(out/'intervals.tsv', mask.intervals, ['accession','start','end','event_id','evidence'])
        write_fasta(out/'reference.fasta', mask.render(reference, settings.get('mode','hard')))


def propose(settings, reference, atlas, output):
    with stage(output, settings, [reference, atlas]) as (out, _):
        policy = settings['policy']
        digest = sha(reference)
        threshold = int(settings.get('min_competing_taxa', 2))
        if threshold < 2:
            raise ValueError('At least two competing taxa are required')
        evidence = rows(atlas)
        for row in evidence:
            if row['rank'] != policy['rank'] or row['reference_sha256'] != digest:
                raise ValueError('Atlas rank/reference differs from mask policy')
        records = [dict(accession=r['accession'], start=r['start'], end=r['end'],
                        event_id=r['evidence_id'], evidence=r.get('method','shared_exact_kmer'))
                   for r in evidence if int(r['taxa_count'])>=threshold and r['unknown_taxonomy']=='false'
                   and (settings.get('include_low_complexity', False) or r['low_complexity']=='false')]
        mask = MaskSet.create(reference, records, policy)
        mask.save(out/'mask.json')
        # A matched negative control preserves each accession's union interval lengths.
        rng = random.Random(int(settings.get('seed',17)))
        control = []
        for accession, length in mask.lengths.items():
            widths = [b-a for a,b in mask.merged(accession)]
            rng.shuffle(widths)
            remaining = length-sum(widths)-max(0,len(widths)-1)
            # Uniform weak composition of unmasked space: never truncate/overlap masks.
            bars = sorted(rng.sample(range(remaining+len(widths)), len(widths))) if widths else []
            cuts = [-1, *bars, remaining+len(widths)]
            gaps = [cuts[i+1]-cuts[i]-1 for i in range(len(cuts)-1)]
            position = 0
            for i,width in enumerate(widths):
                position += gaps[i] + int(i>0)
                control.append(dict(accession=accession,start=position,end=position+width,evidence='random_control'))
                position += width
        MaskSet.create(reference, control, {**policy, 'control':'matched_random'}).save(out/'random_mask.json')
