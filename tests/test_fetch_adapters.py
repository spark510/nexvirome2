import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('fetch_adapters', ROOT/'scripts/fetch_adapters.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture
def document():
    entries = json.loads((ROOT/'configs/adapters/catalog.json').read_text())['entries']
    sequences = {e['id']:e['sequence'] for e in entries}
    return '\n\n'.join([
        'TruSeq single index kits:',
        '* Read 1: '+sequences['truseq_r1']+'\n* Read 2: '+sequences['truseq_r2'],
        'AmpliSeq for Illumina; Nextera XT:', '* '+sequences['nextera_dna_prep'],
        'Illumina Stranded mRNA Prep, Ligation:', '* '+sequences['stranded_rna_ligation'],
        'TruSeq Small RNA:', '* '+sequences['truseq_small_rna']])


def test_offline_export_provenance_and_hashes(tmp_path, document):
    source = tmp_path/'source.md'
    source.write_text(document)
    out = tmp_path/'panel'
    manifest = module.fetch(out, source)
    assert manifest['source_mode'] == 'local_snapshot_unverified'
    for name, digest in manifest['files'].items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest() == digest
    for mate in ('R1','R2'):
        assert (out/f'screen_{mate}.fasta').read_text() == (ROOT/f'configs/adapters/screen_{mate}.fasta').read_text()
    assert all(e['retrieved_date'] is None for e in json.loads((out/'catalog.json').read_text())['entries'])
    with pytest.raises(FileExistsError):
        module.fetch(out, source)


@pytest.mark.parametrize('old,new', [('Read 1:', 'Read 2:'), ('TruSeq Small RNA:', 'Unknown kit:'), ('* Read 1:', '* Read one:')])
def test_layout_changes_fail_before_creating_output(tmp_path, document, old, new):
    source = tmp_path/'source.md'
    source.write_text(document.replace(old,new))
    with pytest.raises(ValueError):
        module.fetch(tmp_path/'panel', source)
    assert not (tmp_path/'panel').exists()
