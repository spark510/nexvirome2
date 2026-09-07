"""Aggregate results without claiming success when baseline errors are absent."""
import json
from pathlib import Path
from .common import dump, stage


def report(settings, run_root, output):
    metric_files = sorted(Path(run_root).glob('*/evaluation/*/metrics.json'))
    if not metric_files:
        raise ValueError('No sample/evaluation/method/metrics.json found')
    with stage(output, settings, metric_files) as (out, manifest):
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        records, source_records, resources = [], [], []
        for path in metric_files:
            for cutoff, m in json.loads(path.read_text()).items():
                records.append({'sample': path.parents[2].name, 'method': path.parent.name, 'min_length': int(cutoff),
                                **{k: v for k, v in m.items() if k != 'per_source'}})
                source_records.extend({'sample': path.parents[2].name, 'method': path.parent.name,
                                       'min_length': int(cutoff), 'accession': a, **values}
                                      for a, values in m.get('per_source', {}).items())
        for path in Path(run_root).glob('*/assemblies/*/manifest.json'):
            data = json.loads(path.read_text())
            resources.append({'sample': path.parents[2].name, 'method': path.parent.name,
                              'elapsed_seconds': data.get('elapsed_seconds'),
                              'max_rss_kb': max((c.get('max_rss_kb',0) for c in data.get('commands',[])), default=0) or None})
        for path in Path(run_root).glob('*/correction/manifest.json'):
            data = json.loads(path.read_text())
            summary = json.loads((path.parent/'summary.json').read_text())
            resources.append({'sample': path.parents[1].name, 'method':'corrected',
                              'elapsed_seconds':data.get('elapsed_seconds'), 'split_count':summary['split_count'],
                              'max_rss_kb':max((c.get('max_rss_kb',0) for c in data.get('commands',[])),default=0) or None})
        pd.DataFrame(source_records).to_csv(out/'source_recovery.tsv',sep='\t',index=False)
        pd.DataFrame(resources).to_csv(out/'resources.tsv',sep='\t',index=False)
        frame = pd.DataFrame(records)
        frame.to_csv(out / 'metrics.tsv', sep='\t', index=False)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        subset = frame[frame.min_length == 0]
        grouped = subset.groupby('method')[['wrong_connections', 'genome_recovery']].mean()
        grouped.wrong_connections.plot.bar(ax=axes[0], title='Mean confirmed wrong connections')
        grouped.genome_recovery.plot.bar(ax=axes[1], title='Mean genome recovery')
        for ax in axes:
            ax.tick_params(axis='x', rotation=30)
        fig.tight_layout()
        fig.savefig(out / 'benchmark.png', dpi=180)
        fig.savefig(out / 'benchmark.svg')
        plt.close(fig)
        checks = []
        for sample in sorted(set(subset['sample'])):
            by_method = subset[subset['sample'] == sample].set_index('method')
            if 'metaspades' not in by_method.index or 'corrected' not in by_method.index:
                continue
            before, after = by_method.loc['metaspades'], by_method.loc['corrected']
            reduction = float(1-after.wrong_connections/before.wrong_connections) if before.wrong_connections else None
            loss = float(before.genome_recovery-after.genome_recovery)
            checks.append({'sample': sample, 'connection_reduction': reduction, 'recovery_loss': loss,
                           'target_met': bool(reduction is not None and reduction >= .5 and loss <= .02)})
        dump(out / 'targets.json', {'comparisons': checks, 'interpretation': 'Per-sample targets; independent panel required for a generalization claim.'})
        (out / 'README.md').write_text('NexVirome2 benchmark\n\nZero baseline errors have undefined reduction. '
            'Use single-source controls and a frozen independent panel before interpreting performance.\n', encoding='utf-8')
