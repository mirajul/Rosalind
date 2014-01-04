#!/usr/bin/env python
'''
A solution to a ROSALIND bioinformatics problem.

Problem Title: Global Alignment with Scoring Matrix and Affine Gap Penalty
Rosalind ID: GAFF
Rosalind #: 096
URL: http://rosalind.info/problems/gaff/
'''


def global_alignment_affine_gap_penalty(v, w, scoring_matrix, sigma, epsilon):
    '''Returns the global alignment score of v and w with constant gap peantaly sigma subject to the scoring_matrix.'''
    from numpy import zeros

    # Initialize the matrices.
    S_lower = zeros((len(v)+1, len(w)+1), dtype=int)
    S_middle = zeros((len(v)+1, len(w)+1), dtype=int)
    S_upper = zeros((len(v)+1, len(w)+1), dtype=int)
    backtrack = zeros((len(v)+1, len(w)+1), dtype=int)

    # Initialize the edges with the given penalties.
    for i in xrange(1, len(v)+1):
        S_lower[i][0] = -sigma - (i-1)*epsilon
        S_middle[i][0] = -sigma - (i-1)*epsilon
        S_upper[i][0] = -10*sigma
    for j in xrange(1, len(w)+1):
        S_upper[0][j] = -sigma - (j-1)*epsilon
        S_middle[0][j] = -sigma - (j-1)*epsilon
        S_lower[0][j] = -10*sigma

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            S_lower[i][j] = max([S_lower[i-1][j] - epsilon, S_middle[i-1][j] - sigma])
            S_upper[i][j] = max([S_upper[i][j-1] - epsilon, S_middle[i][j-1] - sigma])
            middle_scores = [S_lower[i][j], S_middle[i-1][j-1] + scoring_matrix[v[i-1], w[j-1]], S_upper[i][j]]
            S_middle[i][j] = max(middle_scores)
            backtrack[i][j] = middle_scores.index(S_middle[i][j]) + 1

   # Initialize the values of i,j and get the minimum score.
    i,j = len(v), len(w)
    max_score = S_middle[i][j]
    v_aligned, w_aligned = v, w

    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Backtrack to the edge of the matrix starting bottom right.
    while i*j != 0:
        if backtrack[i][j] == 1:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i][j] == 3:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        else:
            i -= 1
            j -= 1

    # Prepend the necessary preceeding indels to get to (0,0).
    for repeat in xrange(i):
        w_aligned = insert_indel(w_aligned, 0)
    for repeat in xrange(j):
        v_aligned = insert_indel(v_aligned, 0)

    return max_score, v_aligned, w_aligned

if __name__ == '__main__':
    from scripts import BLOSUM62, ReadFASTA

    # Parse the two input protein strings.
    s, t = [fasta[1] for fasta in ReadFASTA('data/rosalind_gaff.txt')]

    # Get the alignment score.
    score = map(str, global_alignment_affine_gap_penalty(s, t, BLOSUM62(), 11, 1))

    # Print and save the answer.
    print '\n'.join(score)
    with open('output/096_GAFF.txt', 'w') as output_data:
        output_data.write('\n'.join(score))
