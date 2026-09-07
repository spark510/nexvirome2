import subprocess
import time
from pathlib import Path

class CommandRunner:
    """Injectable process execution with the existing log and resource schema."""
    def __init__(self, execute=None):
        self.execute = execute or subprocess.run

    def run(self, args, output, manifest, name, stdout=None):
        output = Path(output)
        args = [str(x) for x in args]
        entry = {'argv': args, 'name': name}
        manifest['commands'].append(entry)
        start = time.monotonic()
        timed = Path('/usr/bin/time').exists()
        executed = ['/usr/bin/time', '-v', '-o', str(output / f'{name}.resources.txt'), *args] if timed else args
        with open(output / f'{name}.stderr.log', 'w') as err, open(stdout or output / f'{name}.stdout.log', 'w') as out:
            result = self.execute(executed, stdout=out, stderr=err, check=False)
        entry.update(returncode=result.returncode, elapsed_seconds=time.monotonic()-start)
        if timed:
            for line in (output / f'{name}.resources.txt').read_text().splitlines():
                if 'Maximum resident set size' in line:
                    entry['max_rss_kb'] = int(line.rsplit(':', 1)[1])
        if result.returncode:
            raise RuntimeError(f'{name} exited {result.returncode}; see {output / (name + ".stderr.log")}')


    __call__ = run

    def version(self, executable, flag, output, manifest):
        return self.run([executable, flag], output, manifest, Path(executable).name + "_version")
