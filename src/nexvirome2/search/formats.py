from ..domain import SearchHit

HIT_FIELDS = ["query","target","fident","qcov","tcov","evalue","bits"]

def read_search(path):
    result = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip().split('\t')
            if values == HIT_FIELDS:
                continue
            if len(values) != len(HIT_FIELDS):
                raise ValueError('Expected search columns: '+','.join(HIT_FIELDS))
            r = dict(zip(HIT_FIELDS,values))
            for key in HIT_FIELDS[2:]:
                r[key] = float(r[key])
            r = SearchHit(**r).to_record()
            result.append(r)
    return result


