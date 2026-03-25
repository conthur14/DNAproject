
def interface():
    # make the interface which asks for a file_name for the query
    print("-----/ DNA Sequence Algorithm /-----")
    print("This program takes in a file name containg \na query and finds the most similar sequence \nof DNA in the file DNA_sequences.txt")
    print("------------------------------------")
    return input("Name of the file containg query: ")


def computeSim(s, t):
    # technique used Needleman-Wunsch Algorithm
    m = len(s)
    n = len(t)

    # create dynamic programming table, inital state = all zeros
    dp_table = []
    for i in range(m+1):
        row = []
        for j in range(n+1):
            row.append(0)
        dp_table.append(row)

    # first row and col = gap pentalities
    for i in range(1, m+1):
        dp_table[i][0] = dp_table[i-1][0] - 1
    for j in range(1, n+1):
        dp_table[i][0] = dp_table[0][j-1] - 1
    
    #  fill table, if the cell above and the the left are the same, increment the cell by 1, else decrement
    for i in range(1, m + 1):
        for j in range(1, n+1):
            if (s[i-1] == t[j-1]):
                match = dp_table[i-1][j-1] + 1
            else:
                match = dp_table[i-1][j-1] - 1
            delete = dp_table[i-1][j] - 1
            insert = dp_table[i][j-1] - 1

            dp_table[i][j] = max(match, delete, insert)
    
    # number of matching nucleotides
    return dp_table[m][n]


def FindMostSimilarSubSequence(D, t):
    # best similarity
    best_sim = float('-inf')
    # best sequence
    best_seq = None
    
    # iterate through each sequennce in D and find the longest similarlity
    for i in range(len(D)):
        sim = computeSim(D[i], t)
        if sim > best_sim:
            best_sim = sim
            best_seq = D[i]
    
    return best_seq, best_sim


def main():
    # open file_name
    file_name = interface()
    
    # Open file_name and read to f_query
    f_query = open(file_name, "r")
    query = f_query.read()

    # open DNA_sequences.txt and read into sequences
    f_sequences = open("DNA_sequences.txt", "r")
    
    dna_sequences = []
    current_seq = []
    
    # go and parse each line, if it starts with > then it is a name and we can exclude it
    for line in f_sequences:
        line = line.strip()

        if line.startswith(">"):
            if current_seq:
                dna_sequences.append("".join(current_seq))
                current_seq = []
        else: 
            current_seq.append(line)

    if current_seq:
        dna_sequences.append("".join(current_seq))

    #print(dna_sequences)

    # returns the most similar subsequence in dna_sequence to the query
    seq, sim = FindMostSimilarSubSequence(dna_sequences, query)
    print("DNA Query:\n" + str(query) + "\n")
    print("Best Sequence:\n" + str(seq) + "\n")
    print("Number of Similar Nucleotides:\n" + str(sim) + "\n")
    
if __name__ == "__main__":
    main()