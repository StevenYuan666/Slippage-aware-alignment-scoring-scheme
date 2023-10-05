def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = file.read().split('>')[1:]
        sequences = [seq.split('\n', 1)[1].replace('\n', '') for seq in sequences]
    return sequences[0], sequences[1]


def substitution_cost(a, b, match_score, mismatch_score):
    return match_score if a == b else mismatch_score


def slippage_aware_alignment(file_path, match_score, mismatch_score, cs, cn):
    S, T = read_fasta(file_path)
    m, n = len(S), len(T)
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    Trace = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = i * cn
        Trace[i][0] = "UP"
    for j in range(1, n + 1):
        D[0][j] = j * cn
        Trace[0][j] = "LEFT"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            matchScore = D[i - 1][j - 1] + substitution_cost(S[i - 1], T[j - 1], match_score, mismatch_score)
            delScore = D[i - 1][j] + (cs if S[i - 1] == S[i - 2] else cn)
            insScore = D[i][j - 1] + (cs if T[j - 1] == T[j - 2] else cn)

            max_score = max(matchScore, delScore, insScore)
            D[i][j] = max_score
            if max_score == matchScore:
                Trace[i][j] = "DIAG"
            elif max_score == delScore:
                Trace[i][j] = "UP"
            else:
                Trace[i][j] = "LEFT"

    alignmentS, alignmentT = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if Trace[i][j] == "DIAG":
            alignmentS = S[i - 1] + alignmentS
            alignmentT = T[j - 1] + alignmentT
            i -= 1
            j -= 1
        elif Trace[i][j] == "UP":
            alignmentS = S[i - 1] + alignmentS
            alignmentT = "-" + alignmentT
            i -= 1
        else:
            alignmentS = "-" + alignmentS
            alignmentT = T[j - 1] + alignmentT
            j -= 1

    print(f'Optimal Alignment Score: {D[m][n]}')
    print(f'Optimal Alignment:\n{alignmentS}\n{alignmentT}')


if __name__ == '__main__':
    # small
    print("Align Small")
    file_path = 'hw1_brca1_3utr_small.fa'
    match_score = 1
    mismatch_score = -1
    cs = -1
    cn = -2
    slippage_aware_alignment(file_path, match_score, mismatch_score, cs, cn)
    # full
    print()
    print("Align Full")
    file_path = 'hw1_brca1_3utr_full.fa'
    slippage_aware_alignment(file_path, match_score, mismatch_score, cs, cn)
