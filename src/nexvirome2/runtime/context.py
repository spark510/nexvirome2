import subprocess
import time
import json
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from ..io.formats import dump, sha

class RunContext:
    """Own one fresh execution directory and its provenance lifecycle."""
    def __init__(self, output, settings, inputs=()):
        self.output = Path(output)
        self.settings = dict(settings)
        self.inputs = tuple(inputs)
        self.manifest = {}

    def __enter__(self):
        self._manager = self._lifecycle()
        return self._manager.__enter__()

    def __exit__(self, *exc):
        return self._manager.__exit__(*exc)

    @contextmanager
    def _lifecycle(self):
        """Accept a scheduler-created empty directory, never reuse stage artifacts."""
        output, settings, inputs = self.output, self.settings, self.inputs
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        if any(output.iterdir()):
            raise FileExistsError(f'Output directory contains existing artifacts: {output}')
        started = time.monotonic()
        manifest = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
                    'code_revision': None, 'inputs': {}, 'commands': []}
        self.manifest = manifest
        # Exclusive ownership also prevents two workers claiming an empty directory.
        with open(output / 'manifest.json', 'x', encoding='utf-8') as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write('\n')
        try:
            dump(output / 'config.json', settings)
            try:
                manifest['code_revision'] = subprocess.check_output(
                    ['git', 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL, text=True).strip()
            except (OSError, subprocess.CalledProcessError):
                pass
            code_root = Path(__file__).parents[1]
            manifest['code_checksums'] = {str(p.relative_to(code_root)): sha(p) for p in code_root.rglob('*.py')}
            for path in inputs:
                manifest['inputs'][str(path)] = sha(path)
            dump(output / 'manifest.json', manifest)
            yield self
        except BaseException as exc:
            manifest.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            raise
        else:
            manifest['status'] = 'complete'
        finally:
            manifest['elapsed_seconds'] = time.monotonic() - started
            manifest['outputs'] = {str(p.relative_to(output)): sha(p) for p in output.rglob('*')
                                   if p.is_file() and p.name != 'manifest.json'}
            dump(output / 'manifest.json', manifest)


