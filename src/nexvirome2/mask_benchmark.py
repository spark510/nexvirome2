"""Same unmasked-origin queries across all reference masking arms."""
from pathlib import Path
from .common import stage, config, rows, table, dump, sha
from .classification import classify, evaluate
from .decisions import paired_interval


def run(settings, queries, reference, taxonomy, truth, arms_file, output):
    arms = config(arms_file)
    if not isinstance(arms,dict) or 'unmasked' not in arms or arms['unmasked'] is not None:
        raise ValueError('Arms JSON requires an unmasked: null baseline')
    for name in arms:
        if not name or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in name):
            raise ValueError('Unsafe arm name')
    masks = {name:str((Path(arms_file).resolve().parent/path).resolve()) if path else None for name,path in arms.items()}
    with stage(output,settings,[queries,reference,taxonomy,truth,arms_file,*[p for p in masks.values() if p]]) as (out,_):
        units, summaries = {}, []
        for arm,mask in masks.items():
            folder = out/arm
            classify(settings,queries,reference,taxonomy,folder/'classification',mask)
            evaluate(settings,folder/'classification/predictions.tsv',truth,folder/'evaluation',
                     folder/'classification/classification.json')
            units[arm] = {r['unit']:r for r in rows(folder/'evaluation/units.tsv')}
            summaries.append(dict(arm=arm,**config(folder/'evaluation/metrics.json')))
        comparisons = []
        for arm,values in units.items():
            if arm=='unmasked':
                continue
            base = units['unmasked']
            comparisons.append(dict(arm=arm,independent_units=len(values),
                wrong_reduction=paired_interval([float(base[u]['wrong'])-float(values[u]['wrong']) for u in base]),
                correct_change=paired_interval([float(values[u]['correct'])-float(base[u]['correct']) for u in base]),
                unclassified_change=paired_interval([float(values[u]['unclassified'])-float(base[u]['unclassified']) for u in base])))
        table(out/'summary.tsv',summaries,['arm','total','correct','wrong','ambiguous','unclassified','precision','independent_units'])
        dump(out/'comparison.json',dict(comparisons=comparisons,queries_sha256=sha(queries),
             reference_sha256=sha(reference),backend='diagnostic_exact_seed',
             biological_performance_validated=False,
             limitation='Legacy classifier, sequencing-error profiles and real heldout panels require separate validation'))
