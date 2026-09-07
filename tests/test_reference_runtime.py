import json
import pytest
from nexvirome2.common import stage, write_fasta
from nexvirome2.reference import accessions, parse_genbank
from nexvirome2.panel import kmers, low_complexity
from nexvirome2.cli import parser


def test_stage_failure_and_no_overwrite(tmp_path):
    out = tmp_path/'failed'
    with pytest.raises(ValueError):
        with stage(out, {'seed':42}):
            raise ValueError('test failure')
    assert json.loads((out/'manifest.json').read_text())['status'] == 'failed'
    with pytest.raises(FileExistsError):
        with stage(out, {}):
            pass


def test_accessions_must_be_versioned(tmp_path):
    path = tmp_path/'accessions.txt'
    path.write_text('NC_001802.1\nNC_001802.1\n')
    assert accessions(path) == ['NC_001802.1']
    path.write_text('NC_001802\n')
    with pytest.raises(ValueError):
        accessions(path)


def test_genbank_preserves_segment_and_short_sequence(tmp_path):
    path = tmp_path/'x.xml'
    path.write_text('<GBSet><GBSeq><GBSeq_accession-version>NC_000001.1</GBSeq_accession-version>'
                   '<GBSeq_sequence>acgt</GBSeq_sequence><GBSeq_topology>circular</GBSeq_topology>'
                   '<GBSeq_feature-table><GBFeature><GBFeature_key>source</GBFeature_key>'
                   '<GBFeature_quals><GBQualifier><GBQualifier_name>segment</GBQualifier_name>'
                   '<GBQualifier_value>RNA1</GBQualifier_value></GBQualifier></GBFeature_quals>'
                   '</GBFeature></GBSeq_feature-table></GBSeq></GBSet>')
    seq, metadata = parse_genbank(path)
    assert seq['NC_000001.1'] == 'ACGT'
    assert metadata[0]['segment'] == 'RNA1'
    assert metadata[0]['completeness'] == 'unknown'


def test_canonical_kmers():
    from nexvirome2.common import revcomp
    seq = 'ACGTACCTGATCG'
    assert set(kmers(seq,5)) == set(kmers(revcomp(seq),5))
    assert low_complexity('A'*31)


def test_correction_cli_has_no_truth_inputs():
    args = ['correct','--config','c','--contigs','a','--graph','g','--paths','p','--r1','r1','--r2','r2','--output','o']
    assert parser().parse_args(args).command == 'correct'
    with pytest.raises(SystemExit):
        parser().parse_args(args+['--truth','secret'])
