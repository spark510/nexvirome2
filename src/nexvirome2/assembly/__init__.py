"""Keep complete native assembler output, including intermediate graph evidence."""
from ..common import command, dump, stage, version


def assemble(settings, assembler, r1, r2, output):
    with stage(output, settings, [r1, r2]) as (out, manifest):
        native = out / 'native'
        threads = int(settings.get('threads', 4))
        if assembler == 'metaspades':
            executable = settings.get('spades', 'spades.py')
            version(executable, '--version', out, manifest)
            args = [executable, '--meta', '-1', r1, '-2', r2, '-o', native,
                    '-t', threads, '-m', settings.get('memory_gb', 16)]
            result = {'contigs': 'native/contigs.fasta', 'graph': 'native/assembly_graph_with_scaffolds.gfa',
                      'paths': 'native/contigs.paths', 'before_repeat_resolution': 'native/before_rr.fasta'}
        elif assembler == 'megahit':
            executable = settings.get('megahit', 'megahit')
            version(executable, '--version', out, manifest)
            args = [executable, '-1', r1, '-2', r2, '-o', native, '-t', threads]
            result = {'contigs': 'native/final.contigs.fa'}
        else:
            raise ValueError(f'Unsupported assembler: {assembler}')
        command(args, out, manifest, 'assembly')
        if not (out / result['contigs']).is_file():
            raise ValueError('Assembler did not produce expected contig file')
        dump(out / 'result.json', result)
