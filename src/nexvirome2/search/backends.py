from dataclasses import dataclass
from typing import Protocol
from .formats import HIT_FIELDS


class SearchBackend(Protocol):
    name: str
    def validate(self,model,structure_map): ...
    def arguments(self,query,database,raw,temporary,settings,model=None): ...


def easy_search_arguments(name,query,database,raw,temporary,settings):
    return [name,'easy-search',query,database,raw,temporary,'--format-output',','.join(HIT_FIELDS),
            '--threads',settings.get('threads',4),'-e',settings.get('search_evalue',10),
            '--max-seqs',settings.get('max_hits',100)]


@dataclass(frozen=True)
class MMseqsBackend:
    name: str = 'mmseqs'

    def validate(self,model,structure_map):
        if model or structure_map:
            raise ValueError('MMseqs uses amino acid sequences only')

    def arguments(self,query,database,raw,temporary,settings,model=None):
        return easy_search_arguments(self.name,query,database,raw,temporary,settings)


@dataclass(frozen=True)
class FoldseekBackend:
    name: str = 'foldseek'

    def validate(self,model,structure_map):
        if bool(model)==bool(structure_map):
            raise ValueError('Foldseek needs exactly one of --model or --structure-map')

    def arguments(self,query,database,raw,temporary,settings,model=None):
        args=easy_search_arguments(self.name,query,database,raw,temporary,settings)
        if model:
            args.extend(['--prostt5-model',model])
        return args


def backend_for(name):
    backends={'mmseqs':MMseqsBackend,'foldseek':FoldseekBackend}
    if name not in backends:
        raise ValueError('Unsupported search backend')
    return backends[name]()
