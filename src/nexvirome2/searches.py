"""Compatible search and coordinate-tool entry points."""
from pathlib import Path
from .common import command, dump, sha, stage, version
from .runtime.artifacts import artifact_files
from .search.formats import HIT_FIELDS, read_search
from .search.service import SearchService
from .search.backends import backend_for

def search(settings,proteins,database,output,backend='foldseek',model=None,structure_map=None,cache=None):
    return SearchService(backend_for(backend),command,version).run(settings,proteins,database,output,model,structure_map,cache)

def download_database(settings,name,output):
    if name not in ('ProstT5','PDB','BFVD'):
        raise ValueError('Supported requested databases: ProstT5, PDB, BFVD')
    with stage(output,settings) as (out,manifest):
        version('foldseek','version',out,manifest)
        command(['foldseek','databases',name,out/'db',out/'tmp','--threads',settings.get('threads',4)],out,manifest,'download')
        dump(out/'database.json',{'name':name,'prefix':'db','files':{str(p.relative_to(out)):sha(p) for p in artifact_files(out/'db')}})


def motif(settings,structures,output,index=None,query=None,residues=None):
    if query and not residues:
        raise ValueError('Explicit chain/residue motif required')
    inputs = artifact_files(structures) if structures else artifact_files(index)
    if query:
        inputs.append(Path(query))
    with stage(output,settings,inputs) as (out,manifest):
        version('folddisco','--version',out,manifest)
        if structures:
            index = out/'index'
            command(['folddisco','index','-p',structures,'-i',index,'-t',settings.get('threads',4)],out,manifest,'index')
        if query:
            command(['folddisco','query','-i',index,'-p',query,'-q',residues,'-t',settings.get('threads',4)],
                    out,manifest,'motif',stdout=out/'motif_results.tsv')
        dump(out/'result.json',{'index':str(index),'query':str(query) if query else None,'residues':residues,
                                'status':'raw motif evidence; not a viral identity or chimera call'})


def align_structures(settings,structures,output):
    files = [p for p in artifact_files(structures) if p.suffix.lower() in ('.pdb','.cif','.mmcif')]
    if len(files)<2:
        raise ValueError('FoldMason requires at least two coordinate files')
    with stage(output,settings,files) as (out,manifest):
        version('foldmason','version',out,manifest)
        command(['foldmason','easy-msa',*files,out/'alignment',out/'tmp','--threads',settings.get('threads',4),
                 '--report-mode','1'],out,manifest,'align_structures')


