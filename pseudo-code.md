```
function SlippageAwareAlignment(S, T, M, cs, cn):
    m = length(S)
    n = length(T)
    // Initialize scoring matrix D and traceback matrix Trace
    D = m x n matrix initialized to 0
    Trace = m x n matrix initialized to NULL
    
    // Base cases for gaps at start of sequences
    for i from 1 to m:
        D[i][0] = i * cn  // non-slippage gap penalty for starting gaps
        Trace[i][0] = "UP"
    for j from 1 to n:
        D[0][j] = j * cn  // non-slippage gap penalty for starting gaps
        Trace[0][j] = "LEFT"
    
    // Fill in the matrices
    for i from 1 to m:
        for j from 1 to n:
            matchScore = D[i-1][j-1] + M[S[i], T[j]]  // substitution cost from matrix M
            if S[i] == S[i-1]:
                delScore = D[i-1][j] + cs  // slippage gap penalty for deletions
            else:
                delScore = D[i-1][j] + cn  // non-slippage gap penalty for deletions
            if T[j] == T[j-1]:
                insScore = D[i][j-1] + cs  // slippage gap penalty for insertions
            else:
                insScore = D[i][j-1] + cn  // non-slippage gap penalty for insertions
            
            // Determine the best score and traceback direction
            D[i][j] = max(matchScore, delScore, insScore)
            if D[i][j] == matchScore:
                Trace[i][j] = "DIAG"
            elif D[i][j] == delScore:
                Trace[i][j] = "UP"
            else:
                Trace[i][j] = "LEFT"
    
    // Traceback to find the optimal alignment
    alignmentS = ""
    alignmentT = ""
    i = m
    j = n
    while i > 0 or j > 0:
        if Trace[i][j] == "DIAG":
            alignmentS = S[i] + alignmentS
            alignmentT = T[j] + alignmentT
            i -= 1
            j -= 1
        elif Trace[i][j] == "UP":
            alignmentS = S[i] + alignmentS
            alignmentT = "-" + alignmentT
            i -= 1
        else:
            alignmentS = "-" + alignmentS
            alignmentT = T[j] + alignmentT
            j -= 1
    
    return D[m][n], alignmentS, alignmentT

```

