import glob

fasta_files = glob.glob(('../fasta_data/*'))
seqs = []

# Read & Format
for file_name in fasta_files:
        f = open(file_name, "r")
        flines = f.readlines()
        seq = ""
        for el in flines[1:]:
                seq += el.strip()
        seqs.append(seq) # Place desired length of comparison here.
        f.close()

# Var names from Wikipedia pseudocode
d = -5 # Gap penalty
A = seqs[0] # First sequence to be compared
B = seqs[1] # Second ""
I = range(len(seqs[0])) # To help iterate (Pythonic)
J = range(len(seqs[1])) # ""
F = [[0 for i in seqs[0]] for j in seqs[1]] # Fill a 2D array with zeroes
# Similarity matrix from Wikipedia:
S = \
{'A': {'A': 10, 'G': -1, 'C': -3, 'T': -4},
 'G': {'A': -1, 'G':  7, 'C': -5, 'T': -3},
 'C': {'A': -3, 'G': -5, 'C':  9, 'T':  0},
 'T': {'A': -4, 'G': -3, 'C':  0, 'T':  8}}

# Initialization
for i in I:
        F[i][0] = d * i
for j in J:
        F[0][j] = d * j
        
# Scoring
for i in I[1:]:
        for j in J[1:]:
                Match = F[i-1][j-1] + S[A[i]][B[j]]
                Delete = F[i-1][j] + d
                Insert = F[i][j-1] + d
                F[i][j] = max(Match, Insert, Delete)

# Traceback
AlignmentA = ""
AlignmentB = ""
i = len(seqs[0]) - 1
j = len(seqs[1]) - 1

while (i > 0 and j > 0):
        Score = F[i][j]
        ScoreDiag = F[i - 1][j - 1]
        ScoreUp = F[i][j - 1]
        ScoreLeft = F[i - 1][j]
        if (Score == ScoreDiag + S[A[i]][B[j]]):
                AlignmentA = A[i] + AlignmentA
                AlignmentB = B[j] + AlignmentB
                i -= 1
                j -= 1
        elif (Score == ScoreLeft + d):
                AlignmentA = A[i] + AlignmentA
                AlignmentB = "-" + AlignmentB
                i -= 1
        elif (Score == ScoreUp + d):
                AlignmentA = "-" + AlignmentA
                AlignmentB = B[j] + AlignmentB
                j -= 1
        else:
                print("algorithm error?")
while (i > 0):
        AlignmentA = A[i] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i -= 1
while (j > 0):
        AlignmentA = "-" + AlignmentA
        AlignmentB = B[j] + AlignmentB
        j -= 1

# Similarity
lenA = len(AlignmentA)
lenB = len(AlignmentB)
sim1 = ""
sim2 = ""
len0 = 0
k = 0
total = 0
similarity = 0.0

if (lenA > lenB):
        sim1 = AlignmentA
        sim2 = AlignmentB
        len0 = lenA
else:
        sim1 = AlignmentB
        sim2 = AlignmentA
        len0 = lenB

while (k < len0):
        if (sim1[k] == sim2[k]):
                total += 1
        k += 1
        
similarity = total / len0 * 100

print("Alignment A:\n" + AlignmentA)
print("Alignment B:\n" + AlignmentB)
print("Similarity is: " + str("%.2f" % similarity) + "%") # Works only on Python 3?
