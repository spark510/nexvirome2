def split_sequences(sequences, cuts):
    output, mapping = {}, []
    for name, sequence in sequences.items():
        points = [0, *sorted(set(cuts.get(name, []))), len(sequence)]
        if any(not 0 < p < len(sequence) for p in points[1:-1]):
            raise ValueError('Invalid split coordinate')
        for i, (start, end) in enumerate(zip(points, points[1:])):
            label = f'{name}__part_{i+1}' if len(points) > 2 else name
            if label in output:
                raise ValueError('Generated FASTA ID collision')
            output[label] = sequence[start:end]
            mapping.append({'fragment': label, 'original': name, 'start': start, 'end': end})
    return output, mapping


