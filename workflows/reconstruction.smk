"""Run with --configfile specifying contigs, graph, paths, r1, r2, output."""
from nexvirome2.common import config as read_config

ROOT = config['output'].rstrip('/')

rule all:
    input:
        ROOT+'/reconstructed.fasta', ROOT+'/graph_coordinate_map.tsv', ROOT+'/summary.json'

rule reconstruct:
    input:
        contigs=config['contigs'], graph=config['graph'], paths=config['paths'],
        r1=config['r1'], r2=config['r2'], settings=config.get('settings','configs/reconstruction.yaml')
    output:
        fasta=ROOT+'/reconstructed.fasta', mapping=ROOT+'/graph_coordinate_map.tsv', summary=ROOT+'/summary.json'
    threads: config.get('threads',4)
    run:
        from nexvirome2.reconstruction import reconstruct
        settings=read_config(input.settings)
        settings['threads']=threads
        reconstruct(settings,input.contigs,input.graph,input.paths,input.r1,input.r2,ROOT)
