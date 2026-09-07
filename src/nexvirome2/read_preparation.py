"""Paired FASTQ validation/QC and optional auditable host-call subtraction."""
from pathlib import Path
import gzip
import math
from .common import stage,openseq,rows,table,dump,config,sha


def records(path):
    with openseq(path) as handle:
        while True:
            header=handle.readline()
            if not header:
                return
            sequence=handle.readline().rstrip('\r\n')
            plus=handle.readline().rstrip('\r\n')
            quality=handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or not sequence or len(sequence)!=len(quality):
                raise ValueError('Expected complete four-line FASTQ records')
            if set(sequence.upper())-set('ACGTRYSWKMBDHVN') or any(not 33<=ord(c)<=126 for c in quality):
                raise ValueError('Invalid FASTQ sequence or Phred+33 quality')
            header_fields=header[1:].split()
            if not header_fields:
                raise ValueError('FASTQ record needs a nonempty identifier')
            identifier=header_fields[0]
            if plus[1:].strip() and plus[1:].split()[0]!=identifier:
                raise ValueError('FASTQ repeated identifier differs from header')
            yield identifier,(header.rstrip('\r\n'),sequence,plus,quality)


def fragment_id(identifier):
    return identifier[:-2] if identifier.endswith(('/1','/2')) else identifier


def prepare(settings,r1,r2,output,host_calls=None,host_metadata=None):
    """Host calls must be independent competitive assignments, not simulator truth.

    TSV: fragment,status; status host_supported/ambiguous/nonhost/unclassified.
    Metadata JSON: r1_sha256,r2_sha256,method,reference_sha256. Host evidence is an
    explicit adapter contract; this function never infers it from read names.
    All removed pairs are preserved in separate paired FASTQs.
    """
    minimum=int(settings.get('min_length',50))
    quality=float(settings.get('min_mean_quality',20))
    max_n=float(settings.get('max_n_fraction',.1))
    if minimum<1 or not math.isfinite(quality) or not 0<=quality<=93 or not 0<=max_n<=1:
        raise ValueError('Invalid read QC thresholds')
    if bool(host_calls)!=bool(host_metadata) or settings.get('remove_host',False) and not host_calls:
        raise ValueError('Host subtraction requires calls and matching metadata')
    with stage(output,settings,[p for p in (r1,r2,host_calls,host_metadata) if p]) as (out,_):
        calls={}
        if host_calls:
            metadata=config(host_metadata)
            if metadata['r1_sha256']!=sha(r1) or metadata['r2_sha256']!=sha(r2):
                raise ValueError('Host calls belong to different reads')
            if not metadata.get('method') or not metadata.get('reference_sha256') or metadata.get('truth_used',False):
                raise ValueError('Host calls need non-truth method/reference provenance')
            for row in rows(host_calls):
                name=row['fragment']
                if name in calls or row['status'] not in ('host_supported','ambiguous','nonhost','unclassified'):
                    raise ValueError('Duplicate/invalid host call')
                calls[name]=row['status']
        from contextlib import ExitStack
        from itertools import zip_longest
        from collections import Counter
        counts=Counter()
        seen=set()
        with ExitStack() as stack:
            outputs={name:[stack.enter_context(gzip.open(out/f'{name}_R{mate}.fastq.gz','wt',newline='\n'))
                            for mate in (1,2)] for name in ('retained','quality_filtered','host_removed')}
            audit=stack.enter_context(open(out/'pairs.tsv','w',encoding='utf-8',newline=''))
            import csv
            writer=csv.DictWriter(audit,fieldnames=['fragment','destination','reason','host_status'],delimiter='\t')
            writer.writeheader()
            for first,second in zip_longest(records(r1),records(r2)):
                if first is None or second is None or fragment_id(first[0])!=fragment_id(second[0]):
                    raise ValueError('FASTQ mates are missing or out of order')
                if first[0].endswith('/2') or second[0].endswith('/1'):
                    raise ValueError('FASTQ mate suffix conflicts with R1/R2 input')
                name=fragment_id(first[0])
                if name in seen:
                    raise ValueError('Duplicate fragment ID')
                seen.add(name)
                reason=[]
                for _,seq,_,qual in (first[1],second[1]):
                    if len(seq)<minimum:
                        reason.append('short')
                    if sum(ord(c)-33 for c in qual)/len(qual)<quality:
                        reason.append('low_quality')
                    if seq.upper().count('N')/len(seq)>max_n:
                        reason.append('excess_N')
                host_status=calls.get(name,'not_assessed')
                destination='quality_filtered' if reason else 'host_removed' if settings.get('remove_host',False) and host_status=='host_supported' else 'retained'
                if destination=='host_removed':
                    reason=['explicit_host_support']
                for handle,read in zip(outputs[destination],(first[1],second[1])):
                    handle.write('\n'.join(read)+'\n')
                counts[destination]+=1
                writer.writerow(dict(fragment=name,destination=destination,reason=';'.join(sorted(set(reason))),host_status=host_status))
        if set(calls)-seen:
            raise ValueError('Host evidence contains unknown fragments')
        dump(out/'summary.json',dict(total_pairs=sum(counts.values()),counts=dict(counts),
             removed_pairs_preserved=True,ambiguity_alone_causes_host_removal=False,
             limitation='Imported host-call accuracy and host/viral/ERV competition require separate validation'))
