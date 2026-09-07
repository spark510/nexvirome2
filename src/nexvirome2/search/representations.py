"""Database-independent Foldseek query representation cache (per query batch)."""
from pathlib import Path
from ..common import config,sha,stage,write_fasta
from ..decisions import fingerprint


def prepare(cache, sequences, model, model_files, tool_version, settings, command):
    model_path=Path(model).resolve()
    model_root=model_path if model_path.is_dir() else model_path.parent
    identity=dict(schema=1,representation='prostt5_3di',version=tool_version,
                  sequences=dict(sorted(sequences.items())),
                  model={str(Path(p).resolve().relative_to(model_root)):sha(p) for p in model_files},
                  threads=int(settings.get('threads',4)))
    key=fingerprint(identity)
    folder=Path(cache)/key
    manifest_file=folder/'manifest.json'
    hit=False
    if manifest_file.exists():
        manifest=config(manifest_file)
        if manifest.get('status')!='complete':
            raise ValueError('Incomplete representation cache; use a fresh cache location')
        for name,digest in manifest['outputs'].items():
            artifact=(folder/name).resolve()
            if not artifact.is_relative_to(folder.resolve()) or not artifact.is_file() or sha(artifact)!=digest:
                raise ValueError('Representation cache checksum mismatch')
        if config(folder/'config.json')!=identity:
            raise ValueError('Representation cache identity mismatch')
        hit=True
    else:
        with stage(folder,identity,model_files) as (out,manifest):
            write_fasta(out/'queries.faa',dict(sorted(sequences.items())))
            command(['foldseek','createdb',out/'queries.faa',out/'query_db','--prostt5-model',model,
                     '--threads',settings.get('threads',4)],out,manifest,'prediction')
            if any(not (out/(prefix+suffix)).is_file() for prefix in ('query_db','query_db_ss')
                   for suffix in ('','.dbtype','.index')):
                raise RuntimeError('Foldseek did not produce amino acid and 3Di databases')
    return folder/'query_db',dict(key=key,cache_hit=hit,manifest_sha256=sha(manifest_file),
                                 granularity='query_batch',database_independent=True)
