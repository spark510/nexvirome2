"""Compatibility exports; new implementations live in io and runtime."""
from contextlib import contextmanager
from pathlib import Path
from .io.formats import config, dump, sha, openseq, fasta, write_fasta, revcomp, table, rows, union_length
from .runtime import RunContext, CommandRunner

@contextmanager
def stage(output, settings, inputs=()):
    with RunContext(output, settings, inputs) as context:
        yield context.output, context.manifest

def command(args, output, manifest, name, stdout=None):
    return CommandRunner().run(args, output, manifest, name, stdout)

def version(executable, flag, output, manifest):
    command([executable, flag], output, manifest, Path(executable).name + '_version')


