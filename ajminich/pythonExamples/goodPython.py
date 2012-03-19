#/usr/bin/python

def get_longest_seq(array_in):
    '''
    Gets the longest consecutive subsequence in the array.
    A lovely dynamic programming algorithm based on CS161 HW6. Operates in O(n^2) time.
    '''

    S = []

    for i in range(len(array_in)):

        # Perform the dynamic programming algorithm on all previous optimal solutions
        solutions = [list(S[i-j]) + [array_in[i]] for j in range(i+1) if array_in[i] > array_in[i-j]]

        # Include the possibility that the current link is part of a new solution
        solutions.append([array_in[i]])

        # Find the longest solution and consider it the 'optimal' solution for the current array element
        maxSolution = max(enumerate(solutions), key = lambda tup: len(tup[1]))
        S.append(maxSolution[1])

    # Find the longest sequence out of all the sequences
    maxElem = max(enumerate(S), key = lambda tup: len(tup[1]))
    longestSequence = maxElem[1]

    print "Longest strictly increasing sequence contains %d elements." % len(longestSequence)
    return longestSequence

if __name__ == '__main__':

    test_sequence = [2,1,4,8,15,10]

    print get_longest_seq(test_sequence)