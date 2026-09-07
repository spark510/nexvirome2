"""Check executable availability and startup after environment installation."""
import shutil
import subprocess

TOOLS = {
    'spades.py': '--version', 'megahit': '--version', 'art_illumina': '-h',
    'blastn': '-version', 'bowtie2': '--version', 'bowtie2-build': '--version',
    'samtools': '--version', 'snakemake': '--version',
}


def check(extended=False):
    results = []
    tools = {**TOOLS,**({'mmseqs':'version','foldseek':'version','folddisco':'--version','foldmason':'version',
                       'genomad':'--version','checkv':'--version'} if extended else {})}
    for tool, flag in tools.items():
        executable = shutil.which(tool)
        entry = {'tool': tool, 'path': executable, 'ok': False}
        if executable is None:
            entry['error'] = 'not found in PATH'
        else:
            try:
                run = subprocess.run([executable, flag], capture_output=True, text=True,
                                     errors='replace', timeout=30, check=False)
                lines = (run.stdout+'\n'+run.stderr).strip().splitlines()
                entry.update(ok=run.returncode == 0, returncode=run.returncode,
                             banner=' | '.join(lines[:3]))
            except (OSError, subprocess.TimeoutExpired) as exc:
                entry['error'] = str(exc)
        results.append(entry)
    return {'ok': all(row['ok'] for row in results), 'tools': results,
            'note': 'Executable startup checks only; run a Linux pilot to validate biological tool integration.'}
