"""Exact sequence dereplication with source lineage; no genome joining."""
from collections import defaultdict
import hashlib
from pathlib import Path
from .common import fasta, revcomp, rows, table, stage, write_fasta, dump, sha


def build(settings, inputs, output):
    entries=rows(inputs)
    sources={}
    for row in entries:
        if not row['source'] or row['source'] in sources:
            raise ValueError('Catalog source IDs must be unique')
        sources[row['source']]=(Path(inputs).resolve().parent/row['fasta']).resolve()
    with stage(output,settings,[inputs,*sources.values()]) as (out,_):
        representatives,lineage={},[]
        for source,path in sorted(sources.items()):
            for contig,sequence in fasta(path).items():
                if not sequence:
                    raise ValueError('Empty catalog sequence')
                reverse=revcomp(sequence)
                canonical=min(sequence,reverse)
                identifier='nv_'+hashlib.sha256(canonical.encode()).hexdigest()
                if identifier in representatives and representatives[identifier]!=canonical:
                    raise ValueError('Catalog sequence hash collision')
                representatives[identifier]=canonical
                lineage.append(dict(catalog_id=identifier,source=source,contig=contig,length=len(sequence),
                                    orientation='+' if canonical==sequence else '-',source_sha256=sha(path)))
        write_fasta(out/'catalog.fasta',dict(sorted(representatives.items())))
        table(out/'lineage.tsv',lineage,['catalog_id','source','contig','length','orientation','source_sha256'])
        dump(out/'summary.json',dict(input_contigs=len(lineage),representatives=len(representatives),
             method='exact_sequence_and_reverse_complement',joined_genomes=False,
             limitations=['No near-identity clustering','No circular-rotation equivalence']))


def split(settings, membership, output):
    """Hold out entire connected groups across cluster, genome/segment and rank keys."""
    with stage(output,settings,[membership]) as (out,_):
        records=rows(membership)
        fields=settings.get('group_fields',['cluster','genome_group'])
        if not fields:
            raise ValueError('At least one biological grouping field is required')
        parents={}
        def find(x):
            parents.setdefault(x,x)
            while parents[x]!=x:
                parents[x]=parents[parents[x]]
                x=parents[x]
            return x
        def unite(a,b):
            a,b=find(a),find(b)
            parents[max(a,b)]=min(a,b)
        seen={}
        accessions=set()
        for row in records:
            acc=row['accession']
            if not acc or acc in accessions:
                raise ValueError('Duplicate/empty split accession')
            accessions.add(acc)
            find(acc)
            for field in fields:
                value=row[field]
                if not value:
                    raise ValueError('Missing biological group membership')
                key=(field,value)
                if key in seen:
                    unite(acc,seen[key])
                seen[key]=acc
        groups=defaultdict(list)
        for acc in accessions:
            groups[find(acc)].append(acc)
        if len(groups)<2:
            raise ValueError('Need at least two independent groups for heldout evaluation')
        import random
        order=sorted(groups)
        random.Random(int(settings.get('seed',17))).shuffle(order)
        fraction=float(settings.get('heldout_fraction',.3))
        if not 0<fraction<1:
            raise ValueError('Heldout fraction must be between 0 and 1')
        count=min(len(order)-1,max(1,round(len(order)*fraction)))
        heldout=set(order[:count])
        result=[dict(accession=acc,unit=group,split='heldout' if group in heldout else 'development')
                for group in sorted(groups) for acc in sorted(groups[group])]
        table(out/'split.tsv',result,['accession','unit','split'])
        dump(out/'summary.json',dict(group_fields=fields,independent_groups=len(groups),
             heldout_groups=count,annotation_database_leakage_audited=False))
