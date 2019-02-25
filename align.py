from Bio.Alphabet import IUPAC
from Bio.SubsMat.MatrixInfo import blosum100
from Bio import pairwise2

blosum100X0 = {k: 0 if 'X' in k else v for k, v in blosum100.items()}

def print_matrix(matrix):
    letters = list(IUPAC.extended_protein.letters)
    header = '   ' + '| '.join(letters)
    lines = [ header ]
    for i in range(len(letters)):
        line = [ ]
        for j in range(i + 1):
            key1 = (letters[i], letters[j])
            key2 = (letters[j], letters[i])
            if key1 in matrix:
                line.append('{:>2d}'.format(matrix[key1]))
            elif key2 in matrix:
                line.append('{:>2d}'.format(matrix[key2]))
            else:
                line.append('{:>2s}'.format('  '))
        lines.append(letters[i] + ' ' + '|'.join(line))
    lines.append(header)
    print('\n'.join(lines))

def _map_same_sequence(seq1, seq2):
    '''
    map the indices of the same sequence before/after alignment
    '''
    j = 0
    m = [ None ] * len(seq1)
    for i in range(len(seq1)):
        if seq1[i] == '-': continue
        while seq2[j] == '-': j += 1
        assert seq1[i] == seq2[j], 'seq1[{}]={} seq2[{}]={}'.format(i, seq1[i], j, seq2[j])
        m[i] = j
        j += 1
    return m

def _map_alignment(seq2a):
    '''
    map the indices of two aligned sequences, seq1a and seq2a. 
    seq1a and seq2a mostly form 1-to-1 correspondance, but if 
    a residue in seq1a maps to - or X in seq2a, find the nearest
    non- - or X residue to the left and right in seq2a.
    '''
    m_left = [ 0 ] * len(seq2a)
    m_right = [ len(seq2a) - 1 ] * len(seq2a)
    for i in range(len(seq2a)):
        if seq2a[i] not in '-X':
            m_left[i] = i
            m_right[i] = i
        else:
            if i > 0: m_left[i] = m_left[i - 1]
            if (i > 0) and (seq2a[i - 1] in '-X'):
                m_right[i] = m_right[i - 1]
            else:
                j = i + 1
                while (j < len(seq2a)) and (seq2a[j] in '-X'): j += 1
                if j < len(seq2a): m_right[i] = j
    return m_left, m_right

def _alignment_symbol(seq1a, seq2a):
    assert len(seq1a) == len(seq2a)
    symbol = [ '-' ] * len(seq1a)
    for i in range(len(seq1a)):
        if (seq1a[i] == '-') or (seq2a[i] == '-'): continue
        if (seq1a[i] == seq2a[i]) or (seq1a[i] == 'X') or (seq2a[i] == 'X'):
            symbol[i] = '|'
        else:
            symbol[i] = ' '
    return symbol


def align(seq1, seq2):
    '''
    mapping between seq1 and seq2 indices by alignment.
    return a dictionary of the following keys:
        seq1a: aligned seq1
        seq2a: aligned seq2
        map_1_2_left:  map seq1 indices to seq2, if mapped to - or X in
                       seq2a, use the nearest non- - or X to the left in seq2
        map_1_2_right: map seq1 indices to seq2, if mapped to - or X in
                       seq2a, use the nearest non- - or X to the right in seq2
        map_2_1_left:  map seq2 indices to seq1, if mapped to - or X in
                       seq1a, use the nearest non- - or X to the left in seq1
        map_2_1_right: map seq2 indices to seq1, if mapped to - or X in
                       seq1a, use the nearest non- - or X to the right in seq1
        map_1_1a: map seq1 indices to seq1a
        map_1a_1: map seq1a indices to seq1
        map_2_2a: map seq2 indices to seq2a
        map_2a_2: map seq2a indices to seq2
        symbol: alignment symbols between seq1a and seq2a. - for gaps, | for matching
                residues (including X), space for unmatched residues
    usage:
        mapping = align(seq1, seq2)
        assert seq1[i] == seq1a[mapping['map_1_1a'][i]]
        assert seq1[i] == seq2[mapping['map_1_2_left'][i]]
    '''
    assert type(seq1) == type(seq2)
    gap_char = '-' if isinstance(seq1, str) else [ '-' ]
    alignment = pairwise2.align.localds(seq1, seq2, blosum100X0, 
        -10, -0.5, gap_char=gap_char, one_alignment_only=True)
    seq1a, seq2a = alignment[0][0:2]
    mapping = {}
    map1 = _map_same_sequence(seq1, seq1a)
    map2 = _map_same_sequence(seq2a, seq2)
    m_left, m_right = _map_alignment(seq2a)
    mapping['map_1_2_left'] = [map2[m_left[map1[i]]] for i in range(len(seq1))]
    mapping['map_1_2_right'] = [map2[m_right[map1[i]]] for i in range(len(seq1))]
    map3 = _map_same_sequence(seq2, seq2a)
    map4 = _map_same_sequence(seq1a, seq1)
    m_left, m_right = _map_alignment(seq1a)
    mapping['map_2_1_left'] = [map4[m_left[map3[i]]] for i in range(len(seq2))]
    mapping['map_2_1_right'] = [map4[m_right[map3[i]]] for i in range(len(seq2))]
    mapping['seq1a'] = seq1a
    mapping['seq2a'] = seq2a
    mapping['map_1_1a'] = map1
    mapping['map_1a_1'] = map4
    mapping['map_2_2a'] = map3
    mapping['map_2a_2'] = map2
    mapping['symbol'] = _alignment_symbol(seq1a, seq2a)
    return mapping
