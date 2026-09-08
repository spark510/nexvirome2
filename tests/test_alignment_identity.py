import json
import pytest
from nexvirome2.correction.support import pair_evidence
from nexvirome2.correction.policy import decision


@pytest.mark.parametrize('cigar,score,expected,ambiguous',[
    ('10M20D40M',0,0,3), ('10M20N40M',0,0,3),
    ('10M20D40M',-1,3,0), ('50M',0,3,0),
])
def test_alignment_identity_is_not_fragment_endpoints(tmp_path,cigar,score,expected,ambiguous):
    lines=[]
    for i in range(3):
        left,right=81+i,241+i
        for flag,operation,value in [(99,'50M',0),(355,cigar,score)]:
            lines.append(f'q{i}\t{flag}\talt\t{left}\t42\t{operation}\t=\t{right}\t210\t*\t*\tAS:i:{value}\n')
        lines.append(f'q{i}\t147\talt\t{right}\t42\t50M\t=\t{left}\t-210\t*\t*\tAS:i:0\n')
    sam=tmp_path/'input.sam'
    sam.write_text(''.join(lines))
    audit=tmp_path/'audit.json'
    support,count=pair_evidence(sam,{'current':(103,200),'alt':(103,200)},audit_file=audit)
    assert support=={'current':0,'alt':expected}
    assert count==ambiguous
    assert decision(0,support['alt'])==('split' if expected else 'unresolved')
    rows=json.loads(audit.read_text())['fragments']
    assert all('mate_placements' in c for r in rows for c in r['candidates'])
    if ambiguous:
        assert all(r['status']=='ambiguous' for r in rows)
