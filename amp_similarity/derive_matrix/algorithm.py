import difflib
import numpy as np
import itertools

def blosum(seq_array:list, xx_matrix:int):

    for ignore_index, seq_1 in enumerate(seq_array):
        for i, seq_2 in enumerate(seq_array, 0):
            if ignore_index == i:
                continue

            similarity = difflib.SequenceMatcher(None, seq_1, seq_2).ratio()
            if similarity >= xx_matrix/100:
                seq_array.pop(i)

    # if only one sequence left after elimination
    if len(seq_array) == 1:
        return None #, None, None
    
    # print(f"Sequences after elimination: {seq_array}")
    array = np.empty(shape=(len(seq_array), len(seq_array[0])), dtype=str)
    # print(f"Array shape: {array.shape}")

    #pack sequences in numpy array
    for i in range(len(array)):
        array[i] = list(seq_array[i])

    # print(f"Array after packing: {array}")
    seq_combinations = list(itertools.permutations(''.join(array.flatten(order='C')),2))
    seq_letters = []
    for x in ''.join(array.flatten(order='C')):
        if not x in seq_letters:
            seq_letters.append(x)
    seq_letters = sorted(seq_letters)
    # print(f"Sequence letters: {seq_letters}")

    counts = {}
    for c in seq_combinations:
        c = ''.join(c)
        counts[c] = 0

    for i in range(len(array[0])):
        for p in itertools.combinations(array[:, i], 2):
            combination = ''.join(p)
            if combination != combination[::-1]:
                counts[combination] += 1
                counts[combination[::-1]] += 1
            else:
                counts[combination] += 1

    c_values = dict(sorted(counts.items()))

    return c_values


def blossum_from_counts(counts):
    c_matrix = np.zeros(shape=(400,), dtype=np.int64)
    for i, val in enumerate(counts.values()):
        c_matrix[i] = val
    c_matrix = c_matrix.reshape((20, 20))

    # compute the q values
    total = np.sum(np.triu(c_matrix)) # np.triu gives us upper triangle of matrix
    q_matrix = c_matrix / total

    p_values = []
    idx = 0
    # compute p values
    for column in q_matrix.T:
        p_values.append(column[idx] + 1/2*(sum(column[np.arange(len(column)) != idx])))
        idx += 1

    log_odds = np.zeros(shape=(c_matrix.shape), dtype=np.float64)

    for i, j in np.ndindex(q_matrix.shape):
        with np.errstate(divide='ignore'):
            if i == j:
                log_odds[i, j] = 2*np.log2((q_matrix[i, i])/(p_values[i]**2))
            else:
                log_odds[i, j] = 2*np.log2((q_matrix[i, j])/(2*p_values[i]*p_values[j]))

    str_log_odds = np.zeros(shape=log_odds.shape, dtype=object)

    for (i, j), value in np.ndenumerate(log_odds.round()):
        if np.isneginf(value):  # Better way to check for -inf
            str_log_odds[i][j] = '-inf'
        else:
            str_log_odds[i][j] = str(int(value))

    return str_log_odds
