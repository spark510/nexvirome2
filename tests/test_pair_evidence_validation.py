import json
import pytest

from nexvirome2.correction.support import pair_evidence
from nexvirome2.io.sam import measure_sam


def pair(name, target='c', offset=0, extra=0, tags='AS:i:0'):
    left, right = 51+offset, 251+offset
    return ''.join(
        f'{name}\t{flag+extra}\t{target}\t{pos}\t42\t50M\t=\t{mate}\t{span}\t*\t*\t{tags}\n'
        for flag, pos, mate, span in [(99,left,right,250),(147,right,left,-250)])


@pytest.mark.parametrize('margin', [0, -1, float('nan'), float('inf'), True, '1', None])
def test_invalid_margin_rejected_before_reading_input(tmp_path, margin):
    with pytest.raises(ValueError, match='score_margin'):
        pair_evidence(tmp_path/'missing.sam', {}, margin)


@pytest.mark.parametrize('tags,accepted', [
    ('XS:i:0', 0), ('AS:f:nan\tXS:i:0', 0),
    ('AS:i:0\tXS:f:nan', 0), ('AS:i:0\tXS:f:inf', 0),
    ('AS:Z:invalid\tXS:i:0', 0), ('AS:i:0\tXS:Z:invalid', 0),
    ('AS:i:0\tXS:i:0', 0), ('AS:i:0\tXS:i:-1', 1),
    ('AS:i:0', 1),
])
def test_uncertain_competitive_scores_excluded(tmp_path, tags, accepted):
    sam = tmp_path/'input.sam'
    sam.write_text('@SQ\tSN:c\tLN:500\n'+pair('q', tags=tags))
    depth, stats = measure_sam(sam, {'c':'A'*500}, [dict(contig='c',start=0,end=500)])
    assert stats['accepted_fragments'] == accepted
    assert stats['aligned_bases'] == 100*accepted
    assert depth == [.2*accepted]


@pytest.mark.parametrize('flag', [512, 1024])
def test_flagged_fragment_cannot_support_any_path(tmp_path, flag):
    sam = tmp_path/'input.sam'
    # Even when only an alternative placement is flagged, reject the whole
    # fragment instead of making the remaining placement appear unique.
    sam.write_text(''.join(pair(f'q{i}','current',i,flag) + pair(f'q{i}','alt',i)
                           for i in range(3)))
    audit = tmp_path/'audit.json'
    counts, ambiguous = pair_evidence(sam, {'current':(100,200),'alt':(100,200)}, audit_file=audit)
    assert counts == {'current':0,'alt':0}
    assert ambiguous == 0
    assert all(row['status']=='excluded_qc_or_duplicate'
               for row in json.loads(audit.read_text())['fragments'])


def test_interleaved_fragments_preserve_competition_and_endpoint_deduplication(tmp_path):
    records = (pair('tie','current') + pair('unique','alt',10)
               + pair('duplicate','alt',10) + pair('tie','alt')).splitlines(True)
    sam = tmp_path/'input.sam'
    audit = tmp_path/'audit.json'
    # Separate mates as well as alternative placements of the same fragment.
    sam.write_text(''.join(records[::2]+records[1::2]))
    counts, ambiguous = pair_evidence(sam, {'current':(100,200),'alt':(100,200)}, audit_file=audit)
    assert counts == {'current':0,'alt':1}
    assert ambiguous == 1
    rows = {row['fragment']:row for row in json.loads(audit.read_text())['fragments']}
    assert rows['tie']['status'] == 'ambiguous'
    assert {c['path'] for c in rows['tie']['candidates']} == {'current','alt'}


def test_empty_sam_writes_valid_audit(tmp_path):
    sam = tmp_path/'empty.sam'
    sam.write_text('')
    audit = tmp_path/'audit.json'
    assert pair_evidence(sam, {'p':(100,200)}, audit_file=audit) == ({'p':0},0)
    assert json.loads(audit.read_text())['fragments'] == []
