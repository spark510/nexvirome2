"""Independent geNomad/CheckV outputs; never treated as cross-source chimera truth."""
from .common import stage,command,version,dump,sha
from .runtime.artifacts import artifact_files


def run(settings,tool,contigs,database,output):
    if tool not in ('genomad','checkv'):
        raise ValueError('Unsupported quality tool')
    artifacts = artifact_files(database)
    with stage(output,settings,[contigs,*artifacts]) as (out,manifest):
        version(tool,'--version',out,manifest)
        native = out/'native'
        args = [tool,'end-to-end','--threads',settings.get('threads',4),contigs,native,database] if tool=='genomad' else [
            tool,'end_to_end',contigs,native,'-d',database,'-t',settings.get('threads',4)]
        command(args,out,manifest,'quality')
        dump(out/'result.json',{'tool':tool,'native':'native','input_sha256':sha(contigs),
             'interpretation':'Independent viral classification/quality evidence; not assembly correctness or universal viral completeness'})


def download(settings,tool,output):
    if tool not in ('genomad','checkv'):
        raise ValueError('Unsupported quality tool')
    with stage(output,settings) as (out,manifest):
        version(tool,'--version',out,manifest)
        (out/'database').mkdir()
        command([tool,'download-database' if tool=='genomad' else 'download_database',out/'database'],out,manifest,'download')
        files = artifact_files(out/'database')
        dump(out/'database.json',{'tool':tool,'files':{str(p.relative_to(out)):sha(p) for p in files},
                                 'note':'Use the actual extracted database directory shown in the download log for qc run'})
