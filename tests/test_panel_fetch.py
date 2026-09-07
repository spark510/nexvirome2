import json
import random
from pathlib import Path
import pytest
from nexvirome2.common import write_fasta, fasta, table


def test_fetch_validate_and_mismatch(tmp_path, monkeypatch):
    import nexvirome2.reference as reference
    requested = tmp_path/'acc.txt'
    requested.write_text('NC_000001.1\n')
    def download(url, target, payload=None):
        Path(target).write_text('<GBSet><GBSeq><GBSeq_accession-version>NC_000001.1</GBSeq_accession-version>'
                               '<GBSeq_sequence>acgt</GBSeq_sequence></GBSeq></GBSet>')
    monkeypatch.setattr(reference,'download',download)
    monkeypatch.setattr(reference.time,'sleep',lambda _:None)
    reference.fetch({}, requested, tmp_path/'snapshot')
    assert reference.validate(tmp_path/'snapshot')['valid']
    (tmp_path/'snapshot/reference.fasta').write_text('>NC_000001.1\nAAAA\n')
    with pytest.raises(ValueError,match='Checksum'):
        reference.validate(tmp_path/'snapshot')
    requested.write_text('NC_000002.1\n')
    with pytest.raises(ValueError,match='mismatch'):
        reference.fetch({}, requested, tmp_path/'mismatch')
    mismatch = json.loads((tmp_path/'mismatch/reconciliation.json').read_text())
    assert mismatch['missing'] == ['NC_000002.1'] and mismatch['extra'] == ['NC_000001.1']


def test_panel_selects_shared_regions_and_disjoint_control(tmp_path, monkeypatch):
    import nexvirome2.panel as panel
    from nexvirome2.alignment import FIELDS
    rng = random.Random(123)
    seq = lambda n: ''.join(rng.choices('ACGT',k=n))
    shared1,shared2 = seq(250),seq(220)
    references = {'A':seq(400)+shared1+seq(400),'B':seq(400)+shared1+seq(400),
                  'C':seq(400)+shared2+seq(400),'D':seq(400)+shared2+seq(400),
                  'E':seq(1050),'F':seq(1050)}
    ref, meta = tmp_path/'ref.fa',tmp_path/'meta.tsv'
    write_fasta(ref,references)
    table(meta,[{'accession':a,'topology':'linear'} for a in references],['accession','topology'])
    def align(query, subject, output, manifest, name, threads):
        a,b = next(iter(fasta(query))),next(iter(fasta(subject)))
        shared = shared1 if {a,b}=={'A','B'} else shared2
        h = {'qseqid':a,'sseqid':b,'pident':100,'length':len(shared),
             'qstart':401,'qend':400+len(shared),'sstart':401,'send':400+len(shared),
             'bitscore':500,'qseq':shared,'sseq':shared}
        path = output/f'{name}.tsv'
        path.write_text('\t'.join(str(h[f]) for f in FIELDS.split())+'\n')
        return path
    monkeypatch.setattr(panel,'blast',align)
    panel.build({},ref,meta,tmp_path/'panel')
    pairs = json.loads((tmp_path/'panel/panel.json').read_text())['pairs']
    assert [(p['a'],p['b']) for p in pairs] == [('A','B'),('C','D'),('E','F')]
    assert pairs[-1]['control'] is True
