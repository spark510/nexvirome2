from dataclasses import dataclass
from ..io.formats import rows
from ..domain import SequenceInterval

def _read_matrix(path, window_file):
    import numpy as np
    data = rows(path)
    regions = rows(window_file)
    for r in regions:
        r['start'],r['end'] = int(r['start']),int(r['end'])
        if r['start'] < 0 or r['end'] <= r['start']:
            raise ValueError('Invalid window interval')
    if not data or not regions or len({r['window_id'] for r in data}) != len(data):
        raise ValueError('Empty or duplicate matrix rows')
    mapping = {r['window_id']:r for r in data}
    if len({r['window_id'] for r in regions}) != len(regions) or set(mapping) != {r['window_id'] for r in regions}:
        raise ValueError('Window and coverage IDs differ')
    names = [key for key in data[0] if key != 'window_id']
    matrix = np.array([[float(mapping[r['window_id']][n]) for n in names] for r in regions])
    if not names or not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError('Coverage must be finite and nonnegative')
    return regions,names,matrix


def _aggregate(matrix,names,library_file):
    import numpy as np
    libraries = library_file
    groups = {r['sample']:r['biological_sample'] for r in libraries}
    if len(groups) != len(libraries) or set(groups) != set(names) or any(not g for g in groups.values()):
        raise ValueError('Library metadata must match matrix sample columns exactly')
    biological = sorted(set(groups.values()))
    if len(biological) < 3:
        raise ValueError('At least three biological samples are required; seeds are not independent samples')
    return biological,np.column_stack([matrix[:,[i for i,n in enumerate(names) if groups[n]==g]].mean(axis=1) for g in biological])


@dataclass
class SampleSet:
    records: list

    @classmethod
    def load(cls, path):
        return cls(rows(path))

    def aggregate(self, matrix, names):
        return _aggregate(matrix, names, self.records)


@dataclass
class CoverageMatrix:
    regions: list
    samples: list
    values: object

    def __post_init__(self):
        import numpy as np
        self.values = np.asarray(self.values, dtype=float)
        if self.values.shape != (len(self.regions),len(self.samples)):
            raise ValueError('Coverage shape differs from window/sample metadata')
        if not self.samples or len(set(self.samples)) != len(self.samples):
            raise ValueError('Empty or duplicate matrix sample names')
        if len({r['window_id'] for r in self.regions}) != len(self.regions):
            raise ValueError('Duplicate matrix windows')
        if not np.isfinite(self.values).all() or (self.values < 0).any():
            raise ValueError('Coverage must be finite and nonnegative')
        for region in self.regions:
            SequenceInterval(region['start'],region['end'])

    @classmethod
    def load(cls, path, window_file):
        return cls(*_read_matrix(path, window_file))

    def aggregate(self, sample_set):
        names, values = sample_set.aggregate(self.values, self.samples)
        return CoverageMatrix(self.regions, names, values)


def read_matrix(path, window_file):
    matrix = CoverageMatrix.load(path, window_file)
    return matrix.regions, matrix.samples, matrix.values


def aggregate(matrix, names, library_file):
    return SampleSet.load(library_file).aggregate(matrix, names)
