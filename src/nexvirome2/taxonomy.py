"""Frozen accession-to-rank taxonomy adapter; absent ranks stay unknown."""
from dataclasses import dataclass
from .common import rows, table, stage, dump, sha


@dataclass(frozen=True)
class TaxonomyResolver:
    records: dict

    @classmethod
    def load(cls, path):
        records = {}
        for row in rows(path):
            accession = row['accession']
            if accession in records:
                raise ValueError(f'Duplicate taxonomy accession: {accession}')
            records[accession] = row
        if not records:
            raise ValueError('Empty taxonomy snapshot')
        return cls(records)

    def taxon(self, accession, rank):
        if accession not in self.records:
            raise ValueError(f'Accession missing from taxonomy: {accession}')
        value = self.records[accession].get(rank, '')
        return None if value in ('', 'unknown', 'NA', '0') else value


def build(settings, mapping, nodes, output, merged=None, deleted=None):
    """Read local NCBI taxdump files, never infer ranks from organism names."""
    with stage(output, settings, [p for p in (mapping, nodes, merged, deleted) if p]) as (out, _):
        def records(path):
            with open(path, encoding='utf-8') as handle:
                for line in handle:
                    yield [part.strip() for part in line.split('|')]
        tree = {}
        for row in records(nodes):
            if row[0] in tree:
                raise ValueError('Duplicate taxdump node')
            tree[row[0]] = (row[1], row[2])
        redirects = {r[0]: r[1] for r in records(merged)} if merged else {}
        removed = {r[0] for r in records(deleted)} if deleted else set()
        ranks = settings.get('ranks', ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom'])
        if not ranks or len(set(ranks)) != len(ranks) or set(ranks) & {'accession','taxid','resolved_taxid','status'}:
            raise ValueError('Invalid taxonomy ranks')
        result, accessions = [], set()
        for source in rows(mapping):
            accession, original = source['accession'], source['taxid']
            if accession in accessions:
                raise ValueError('Duplicate accession in taxonomy mapping')
            accessions.add(accession)
            node, seen = original, set()
            while node in redirects:
                if node in seen:
                    raise ValueError('Cycle in merged taxids')
                seen.add(node)
                node = redirects[node]
            resolved = node
            row = dict(accession=accession, taxid=original, resolved_taxid=resolved,
                       status='deleted' if node in removed else 'known' if node in tree else 'unknown')
            lineage = set()
            while node in tree and node not in removed:
                if node in lineage:
                    raise ValueError('Cycle in taxonomy tree')
                lineage.add(node)
                parent, rank = tree[node]
                if rank in ranks:
                    row[rank] = node
                if parent == node:
                    break
                node = parent
            result.append(row)
        table(out/'taxonomy.tsv', result, ['accession','taxid','resolved_taxid','status', *ranks])
        dump(out/'snapshot.json', {'schema': 1, 'records': len(result), 'ranks': ranks,
                                  'nodes_sha256': sha(nodes)})
