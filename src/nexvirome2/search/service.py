import hashlib
import json
import shutil
from pathlib import Path
from ..common import dump, fasta, rows, sha, stage, table, write_fasta
from ..runtime import CommandRunner
from ..runtime.artifacts import artifact_files
from .formats import HIT_FIELDS, read_search
from .cache import SearchCache
from .backends import SearchBackend

class SearchService:
    """Compose a backend, process runner and persistent cache."""
    def __init__(self, backend: SearchBackend, runner=None, version=None, cache_factory=SearchCache):
        self.backend = backend
        runner = runner or CommandRunner()
        self.command = runner
        self.version = version or runner.version
        self.cache_factory = cache_factory

    def run(self,settings,proteins,database,output,model=None,structure_map=None,cache=None):
        backend = self.backend.name
        self.backend.validate(model, structure_map)
        command, version = self.command, self.version
        db_files = artifact_files(database)
        model_files = artifact_files(model) if model else []
        structures = rows(structure_map) if structure_map else []
        mapping = {}
        for r in structures:
            for key in ('protein_id','query_id','path'):
                if not r.get(key):
                    raise ValueError('Structure map needs protein_id, query_id, path')
            if r['protein_id'] in mapping:
                raise ValueError('One structure chain per protein is required')
            r['path'] = str((Path(structure_map).resolve().parent/r['path']).resolve())
            mapping[r['protein_id']] = r
        input_files = [proteins,*db_files,*model_files]+([structure_map] if structure_map else [])+[r['path'] for r in structures]
        with stage(output,settings,input_files) as (out,manifest):
            sequences = fasta(proteins)
            if structure_map and set(sequences) != set(mapping):
                raise ValueError('Structure map must cover exactly the selected proteins; subset proteins first')
            if len({r['query_id'] for r in structures}) != len(structures):
                raise ValueError('Duplicate structure query IDs')
            version(backend,'version',out,manifest)
            tool_version = (out/f'{backend}_version.stdout.log').read_text().strip()
            provenance = {'backend':backend,'version':tool_version,'settings':settings,
                          'database':{str(p):sha(p) for p in db_files},'model':{str(p):sha(p) for p in model_files},
                          'representation':'coordinates' if structures else 'predicted_3di' if model else 'amino_acids'}
            digest = hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
            keys,fresh,found = {},{},[]
            with self.cache_factory(cache or out/'search_cache.sqlite') as db:
                for name,seq in sequences.items():
                    extra = sha(mapping[name]['path'])+'|'+mapping[name]['query_id'] if structures else ''
                    key = hashlib.sha256((digest+'|'+seq+'|'+extra).encode()).hexdigest()
                    keys[name] = key
                    cached = db.get(key)
                    if cached is not None:
                        found.extend({**r,'query':name} for r in cached)
                    else:
                        fresh[name] = seq
                if fresh:
                    query = out/'queries.faa'
                    aliases = {}
                    if structures:
                        query = out/'structures'
                        query.mkdir()
                        filenames = set()
                        for name in fresh:
                            r = mapping[name]
                            basename = Path(r['path']).name
                            if basename in filenames:
                                raise ValueError('Structure files must have unique basenames')
                            filenames.add(basename)
                            shutil.copyfile(r['path'],query/basename)
                            aliases[r['query_id']] = name
                    else:
                        write_fasta(query,fresh)
                    raw = out/'raw_hits.tsv'
                    if backend=='foldseek' and model and settings.get('representation_cache'):
                        from .representations import prepare
                        query_db,identity=prepare(settings['representation_cache'],fresh,model,model_files,
                                                  tool_version,settings,command)
                        manifest['representation_cache']=identity
                        result_db=out/'alignment_db'
                        command(['foldseek','search',query_db,database,result_db,out/'tmp','-a',
                                 '--threads',settings.get('threads',4),'-e',settings.get('search_evalue',10),
                                 '--max-seqs',settings.get('max_hits',100)],out,manifest,'search')
                        command(['foldseek','convertalis',query_db,database,result_db,raw,
                                 '--format-output',','.join(HIT_FIELDS)],out,manifest,'convert_hits')
                    else:
                        args = self.backend.arguments(query,database,raw,out/'tmp',settings,model)
                        command(args,out,manifest,'search')
                    grouped = {q:[] for q in fresh}
                    for r in read_search(raw):
                        query_id = aliases.get(r['query'],r['query'])
                        if query_id not in grouped:
                            raise ValueError('Unknown search query; check structure chain query_id mapping')
                        grouped[query_id].append({**r,'query':query_id})
                    for name,result in grouped.items():
                        db.put(keys[name],result)
                        found.extend(result)
                table(out/'hits.tsv',found,HIT_FIELDS)
                manifest.update(cache_hits=len(sequences)-len(fresh),searched_proteins=len(fresh),**provenance)
                dump(out/'search_identity.json',provenance)


