"""Expand a frozen panel into mixtures and matched single-source experiments."""
import json
from pathlib import Path
from .common import dump, stage


def prepare(settings, panel_file, output):
    panel_file = Path(panel_file).resolve()
    panel = json.loads(panel_file.read_text())
    with stage(output, settings, [panel_file]) as (out, manifest):
        experiments = {}
        partition = settings.get('partition', 'development')
        excluded = set(settings.get('exclude_accessions', []))
        if partition == 'evaluation' and not excluded:
            raise ValueError('Evaluation requires excluded development accessions')
        for index, pair in enumerate(panel['pairs']):
            if {pair['a'], pair['b']} & excluded:
                raise ValueError('Panel overlaps development references')
            for coverage in settings.get('coverages', [50]):
                for ratio in settings.get('ratios', [1]):
                    if coverage <= 0 or ratio <= 0:
                        raise ValueError('Coverage and ratio must be positive')
                    specs = [{'accession': a, 'fasta': str(panel_file.parent / 'sources' / f'{a}.fasta'),
                              'coverage': coverage/(ratio if i else 1), 'seed_offset': i,
                              'topology': panel.get('metadata', {}).get(a, {}).get('topology', 'linear')}
                             for i, a in enumerate((pair['a'], pair['b']))]
                    for seed in settings.get('seeds', [20260905, 20260906, 20260907]):
                        prefix = f'p{index}_c{coverage}_r{ratio}_s{seed}'
                        for suffix, sources in [('mix', specs), ('single_a', specs[:1]), ('single_b', specs[1:])]:
                            name = prefix+'_'+suffix
                            dump(out / f'{name}.sources.json', {'sources': sources})
                            experiments[name] = {'sources': str((out / f'{name}.sources.json').resolve()),
                                                 'seed': seed, 'partition': partition, 'pair': index,
                                                 'coverage': coverage, 'ratio': ratio, 'sample_type': suffix}
        dump(out / 'experiments.json', {'experiments': experiments})
        manifest['experiments'] = len(experiments)
