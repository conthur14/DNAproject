
def interface():
    # make the interface which asks for a file_name for the query
    print("-----/ DNA Sequence Algorithm /-----")
    print("Fill with Details about the algorithm")
    return input("Name of the File: ")

def computeSim(s, t):
    # technique used Needleman-Wunsch Algorithm
    m = len(s)
    n = len(t)

    # create dynamic programming table
    dp_table = []
    for i in range(m+1):
        row = []
        for j in range(n+1):
            row.append(0)
        dp_table.append(row)

    max_len = 0 
    end_index_s = 0



def FindMostSimilarSubSequence(t, D):
    # best similarity
    best_sim = float('-inf')
    # best sequence
    best_seq = None
    
    # iterate through each sequennce in D and find the longest similarlity
    for i in range(len(D)):
        sim = computeSim(D[i], t)
        if sim > best_sim:
            best_sim = sim
            best_seq = best_seq
    
    return best_seq, best_sim

def main():
    # open file_name
    file_name = interface()
    
    # Open file_name and read to f_query
    f_query = open(file_name, "r")
    query = f_query.read()

    # open DNA_sequences.txt and read into sequences
    f_sequences = open("DNA_sequences.txt", "r")
    sequences = f_sequences.read()

    # parse the lines from the DNA_Sequences.txt file into a list
    lines = sequences.splitlines()
    dna_sequences = []
    for line in lines:
        line = line.strip().upper()
        if all(c in "ACTG" for c in line) and (line != ""):
            dna_sequences.append(line)

    # returns the most similar subsequence in dna_sequence to the query
    sim, seq = FindMostSimilarSubSequence(query, dna_sequences)
    print("Best Sim: " + sim)
    print("Best Sequence" + seq)
    
if __name__ == "__main__":
    main()