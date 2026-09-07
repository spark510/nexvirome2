from pathlib import Path

def artifact_files(path):
    path = Path(path).resolve()
    if path.is_dir():
        files = sorted(p for p in path.rglob('*') if p.is_file())
    else:
        files = sorted(p for p in path.parent.glob(path.name+'*') if p.is_file())
    if not files:
        raise ValueError(f'No database/model artifacts: {path}')
    return files


