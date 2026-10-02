def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    pointer = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize the first col. Fill the first col by aligning seq1 prefixes to gaps.
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + scoring_function(seq1[i - 1], "-")
        pointer[i][0] = "up"

    # Fill the first row by aligning seq2 prefixes to gaps.
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + scoring_function("-", seq2[j - 1])
        pointer[0][j] = "left"

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = score[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = score[i][j - 1] + scoring_function("-", seq2[j - 1])

            score[i][j] = max(diagonal, up, left)

            if score[i][j] == diagonal:
                pointer[i][j] = "diag"
            elif score[i][j] == up:
                pointer[i][j] = "up"
            else:
                pointer[i][j] = "left"

    # Trace back from the bottom right cell
    aln1, aln2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        move = pointer[i][j]
        if move == "diag":
            aln1.append(seq1[i - 1])
            aln2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            aln1.append(seq1[i - 1])
            aln2.append("-")
            i -= 1
        else:
            aln1.append("-")
            aln2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aln1)), "".join(reversed(aln2)), float(score[n][m])
    


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
  
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    pointer = [[None] * (m + 1) for _ in range(n + 1)]

    # The first row and column stay at zero

    best = 0
    best_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = score[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = score[i][j - 1] + scoring_function("-", seq2[j - 1])

            score[i][j] = max(diagonal, up, left, 0)

            if score[i][j] == 0:
                pointer[i][j] = None
            elif score[i][j] == diagonal:
                pointer[i][j] = "diag"
            elif score[i][j] == up:
                pointer[i][j] = "up"
            else:
                pointer[i][j] = "left"

            if score[i][j] > best:
                best = score[i][j]
                best_pos = (i, j)

    # Trace back from the best cell until a zero score
    aln1, aln2 = [], []
    i, j = best_pos
    while i > 0 and j > 0 and score[i][j] > 0:
        move = pointer[i][j]
        if move == "diag":
            aln1.append(seq1[i - 1])
            aln2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            aln1.append(seq1[i - 1])
            aln2.append("-")
            i -= 1
        else:
            aln1.append("-")
            aln2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aln1)), "".join(reversed(aln2)), float(best)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
