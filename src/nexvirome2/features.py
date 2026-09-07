"""Feature application facade with compatible CLI entry points."""
import json
from collections import defaultdict
from pathlib import Path
from .common import command, dump, fasta, rows, sha, stage, table, version
from .coverage.math import windows, composition, cosine, merged
from .coverage.matrix import read_matrix, aggregate
from .io.sam import blocks, sam_records, measure_sam

def neighborhoods(settings,contigs,graph_file,paths_file,output):
    from .graph import Graph,contig_paths,candidates
    with stage(output,settings,[contigs,graph_file,paths_file]) as (out,manifest):
        sequences = fasta(contigs)
        proposed = candidates(Graph.load(graph_file),sequences,contig_paths(paths_file),settings.get('max_repeat_nodes',20))
        flank = int(settings.get('neighborhood_bp',2000))
        if flank < 1:
            raise ValueError('neighborhood_bp must be positive')
        intervals = defaultdict(list)
        for candidate in proposed:
            if candidate['usable']:
                cut = candidate['cut']
                intervals[candidate['contig']].append((max(0,cut-flank-candidate['repeat_length']),min(len(sequences[candidate['contig']]),cut+flank)))
        selected = [{'contig':name,'start':a,'end':b} for name,items in intervals.items() for a,b in merged(items)]
        table(out/'regions.tsv',selected,['contig','start','end'])
        dump(out/'candidates.json',proposed)
        manifest['usable_candidates'] = sum(c['usable'] for c in proposed)


def coverage(settings, contigs, samples_file, output, region_file=None):
    sample_rows = rows(samples_file)
    if not sample_rows or any(not s.get('sample') or not s.get('biological_sample') for s in sample_rows):
        raise ValueError('samples.tsv requires sample and biological_sample columns')
    if len({s['sample'] for s in sample_rows}) != len(sample_rows):
        raise ValueError('Duplicate sample IDs')
    base = Path(samples_file).resolve().parent
    inputs = [contigs,samples_file] + ([region_file] if region_file else [])
    for s in sample_rows:
        if s.get('sam') and s.get('target_sha256') != sha(contigs):
            raise ValueError('Imported SAM requires target_sha256 matching the common contig FASTA')
        for key in (('sam',) if s.get('sam') else ('r1','r2')):
            if not s.get(key):
                raise ValueError('Each sample requires sam or both r1 and r2')
            s[key] = str((base/s[key]).resolve())
            inputs.append(s[key])
    with stage(output,settings,inputs) as (out,manifest):
        target = fasta(contigs)
        selected = None
        if region_file:
            selected = defaultdict(list)
            for r in rows(region_file):
                a,b = int(r['start']),int(r['end'])
                if r['contig'] not in target or not 0 <= a < b <= len(target[r['contig']]):
                    raise ValueError('Invalid requested region')
                selected[r['contig']].append((a,b))
        regions = windows(target,settings.get('window_bp',500),settings.get('step_bp',500),selected)
        if not regions:
            raise ValueError('No selected windows')
        table(out/'windows.tsv',regions,['window_id','contig','start','end'])
        comps = []
        for r in regions:
            seq = target[r['contig']][r['start']:r['end']]
            comps.append({'window_id':r['window_id'],'gc':(seq.count('G')+seq.count('C'))/len(seq),
                          'ambiguous_fraction':sum(x not in 'ACGT' for x in seq)/len(seq),
                          'kmers':composition(seq)})
        dump(out/'composition.json',comps)
        values, normalized, libraries = {},{},[]
        if any(not s.get('sam') for s in sample_rows):
            version('bowtie2','--version',out,manifest)
            command(['bowtie2-build',contigs,out/'target_index'],out,manifest,'index')
        version('samtools','--version',out,manifest)
        for index,sample in enumerate(sample_rows):
            source = sample.get('sam')
            if not source:
                source = out/f'sample_{index}.sam'
                command(['bowtie2','-x',out/'target_index','-1',sample['r1'],'-2',sample['r2'],
                         '-k','2','--very-sensitive','--no-mixed','--no-discordant',
                         '-p',settings.get('threads',4),'-S',source],out,manifest,f'map_{index}')
            sorted_sam = out/f'sample_{index}.name_sorted.sam'
            command(['samtools','sort','-n','-O','SAM','-o',sorted_sam,source],out,manifest,f'sort_{index}')
            depth,stat = measure_sam(sorted_sam,target,regions,settings.get('min_mapq',20))
            values[sample['sample']] = depth
            normalized[sample['sample']] = [v*1e6/stat['aligned_bases'] if stat['aligned_bases'] else 0 for v in depth]
            libraries.append({'sample':sample['sample'],'biological_sample':sample['biological_sample'],**stat})
        names = [s['sample'] for s in sample_rows]
        for filename,matrix in [('coverage.tsv',values),('normalized_coverage.tsv',normalized)]:
            table(out/filename,[{'window_id':r['window_id'],**{n:matrix[n][i] for n in names}} for i,r in enumerate(regions)],['window_id',*names])
        table(out/'libraries.tsv',libraries,['sample','biological_sample','accepted_fragments','excluded_fragments','aligned_bases'])
        manifest.update(target_sha256=sha(contigs),normalization='mean depth / accepted aligned bases * 1e6')


def propose(settings,coverage_file,window_file,library_file,composition_file,output):
    with stage(output,settings,[coverage_file,window_file,library_file,composition_file]) as (out,manifest):
        regions,names,matrix = read_matrix(coverage_file,window_file)
        names,matrix = aggregate(matrix,names,library_file)
        comp = {r['window_id']:r['kmers'] for r in json.loads(Path(composition_file).read_text())}
        if set(comp)!={r['window_id'] for r in regions}:
            raise ValueError('Composition windows differ from coverage windows')
        grouped = defaultdict(list)
        for i,r in enumerate(regions):
            grouped[r['contig']].append(i)
        proposals = []
        for contig,indices in grouped.items():
            indices.sort(key=lambda i:regions[i]['start'])
            for a,b in zip(indices,indices[1:]):
                left,right = regions[a],regions[b]
                if right['start']>left['end']:
                    continue
                cov = cosine(matrix[a],matrix[b])
                ca,cb = comp[left['window_id']],comp[right['window_id']]
                words = sorted(set(ca)|set(cb))
                similarity = cosine([ca.get(w,0) for w in words],[cb.get(w,0) for w in words])
                for method,value in [('coverage',cov),('composition',similarity)]:
                    if value is not None and 1-value>=settings.get(method+'_change',0.5):
                        proposals.append({'contig':contig,'start':min(left['start'],right['start']),
                                          'end':max(left['end'],right['end']),'method':method,'score':max(0,min(1,1-value)),
                                          'detail':json.dumps({'left':left['window_id'],'right':right['window_id']})})
        table(out/'proposals.tsv',proposals,['contig','start','end','method','score','detail'])


