class DnaSeq:
    def __init__(self, accession, seq):
        if not accession or not seq:
            raise ValueError("Accession and sequence can not be empty")
        self.accession = accession
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f"<DnaSeq accession={self.accession}>"


def read_dna(filename):
    sequences = []
    current_accession = None
    current_seq = []
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            if line.startswith('>'):
                if current_accession is not None:
                    # Only add if we have both accesssion and sequence 
                    if current_seq:
                        sequences.append(DnaSeq(current_accession, ''.join(current_seq)))
                current_accession = line[1:]  # Remove '>' character
                current_seq = []
            else:
                current_seq.append(line.upper())  # Store sequence in uppercase

    # Add the last sequence if it exists
    if current_accession is not None and current_seq:
        sequences.append(DnaSeq(current_accession, ''.join(current_seq)))
    
    return sequences


def check_exact_overlap(a, b, min_length=10):
    """
    Find the longest exact overlap between suffix of a and prefix of b.
    Returns lenght of overlap if found (>= min_lenght) otherwise 0.
    """

    max_possible = min(len(a), len(b))
    for L in range(max_possible, min_length - 1, -1):
        if a.seq[-L:] == b.seq[:L]:
            return L
    return 0    


def overlaps(seq_list, overlap_func):
    result = {}
    for a in seq_list:
        overlaps_for_a = {}
        for b in seq_list:
            if a.accession == b.accession:
                continue      # Don't compare with self
            overlap_len = overlap_func(a, b)
            if overlap_len > 0:
                overlaps_for_a[b.accession] = overlap_len
        if overlaps_for_a: # Only add if there are any overlaps
            result[a.accession] = overlaps_for_a
    return result           


#
# Testing code. You should not change any line after this one!
#
def test_class_DnaSeq():
    s1 = DnaSeq('s1', 'ACGT')
    s2 = DnaSeq('s2', 'ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAAT')
    assert len(s1) == 4, 'Your length method (__len__) is not correct.'
    assert len(s2) == 70, 'Your length method (__len__) is not correct.'

    assert str(s1) == '<DnaSeq accession=s1>', 'The __str__ method is not following the specification.'
    assert str(s2) == '<DnaSeq accession=s2>', 'The __str__ method is not following the specification.'

    # The rest of this function is verifying that we are indeed raising an exception.
    status = 0
    try:
        s3 = DnaSeq('', 'ACGT')
    except ValueError:
        status += 1
    try:
        s3 = DnaSeq('s3', None)
    except ValueError:
        status += 1

    try:
        s3 = DnaSeq(None, '')
    except ValueError:
        status += 1
    if status != 3:
        raise Exception('class DnaSeq does not raise a ValueError '
                        'exception with initialised with empty '
                        'accession and sequence.')
    print('DnaSeq passed')


def test_reading():
    dna1 = read_dna('ex1.fa')
    assert len(dna1) == 6, 'The file "ex1.fa" has exactly 6 sequences, but your code does not return that.'
    assert list(map(lambda x: x.accession, dna1)) == [f's{i}' for i in range(6)], 'The accessions are not read correctly'
    print('read_dna passed')

def test_overlap():
    s0 = DnaSeq('s0', 'AAACCC')
    s1 = DnaSeq('s1', 'CCCGGG')
    s2 = DnaSeq('s2', 'TTTTCC')
    data1 = [s0, s1, s2]
    assert check_exact_overlap(s0, s1, 2) == 3
    assert check_exact_overlap(s0, s1) == 0
    assert check_exact_overlap(s1, s2, 2) == 0
    assert check_exact_overlap(s2, s1, 2) == 2

    res0 = overlaps(data1, lambda s1, s2: check_exact_overlap(s1, s2, 2))
    assert len(res0) == 2, 'You get the wrong number of overlaps'
    assert res0 == {'s0': {'s1': 3}, 's2': {'s1': 2}}

    dna_data = read_dna('ex1.fa')
    res1 = overlaps(dna_data, check_exact_overlap)
    assert len(res1) == 5
    for left, right in [('s0', 's1'), ('s1', 's2'), ('s2', 's3'), ('s3', 's4'), ('s4', 's5')]:
        assert res1[left][right], f'Missing overlap of {left} and {right} (in that order)'
    print('overlap code passed')



def test_all():
    test_class_DnaSeq()
    test_reading()
    test_overlap()
    print('Yay, all good')

if __name__ == "__main__":
    test_all()    