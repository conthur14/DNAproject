===========================================================
  DNA Sequence Matching - CSCI 311 Spring 2026 - Team 6
===========================================================

TEAM MEMBERS
------------
  Connor Thurston (cmt012)
  AJ Mehta (ahm057)
  Genesis Evbenaye (oae001)
  Anna Barrios (aeb036)

LANGUAGE
--------
  Python 3 (tested on Python 3.10+)

PREREQUISITES
-------------
  No external libraries are required. This project uses
  only the Python standard library. Any machine with
  Python 3.10 or higher installed can run it.

  To verify your Python version:
      python3 --version

FILES
-----
  README.txt        - This file.
  main.py           - Entry point. Runs the interface and
                      orchestrates the search.
  file_handler.py   - Parses FASTA-formatted sequence files
                      and plain query files.
  algorithms.py     - Contains all similarity/alignment
                      algorithm implementations.
  DNA_sequences.txt - Sample database of known sequences.
  DNA_query.txt     - Sample query sequence.

HOW TO RUN
----------
  1. Open a terminal and navigate to the project directory:

         cd DNAproject

  2. Run the program:

         python3 main.py

  3. The program will prompt you for:
       a) The path to the query sequence file
          (e.g., DNA_query.txt)
       b) The path to the database sequence file
          (e.g., DNA_sequences.txt)

  4. A menu will appear listing the available algorithms:

         =============================================
          DNA Sequence Matching
         =============================================
          1. Longest Common Substring
          2. Longest Common Subsequence (LCS)
          3. Edit Distance
          4. Needleman-Wunsch Alignment
          5. Exit
         =============================================

  5. Select an algorithm by entering its number (1-4).
     The program will compare the query against every
     sequence in the database and report the most
     similar match, including the sequence header
     and the similarity score.

  6. After results are displayed, the menu reappears
     so you can try another algorithm or exit.

EXAMPLE SESSION
---------------
  $ python3 main.py
  Enter path to query sequence file: DNA_query.txt
  Enter path to database sequence file: DNA_sequences.txt

  =============================================
   DNA Sequence Matching
  =============================================
   1. Longest Common Substring
   2. Longest Common Subsequence (LCS)
   3. Edit Distance
   4. Needleman-Wunsch Alignment
   5. Exit
  =============================================
  Select an algorithm (1-5): 2

  Searching with: Longest Common Subsequence (LCS)
  Comparing against 10 sequences...
  [1/10] NC_000011.10 ... done
  [2/10] NT_176377.1 ... done
  ...

  ============ RESULT ============
  Most similar sequence:
  >V01243.1 Rat gene for insulin 2
  Score: 951 (LCS length)
  =================================

  Select an algorithm (1-5): 5
  Goodbye!

INPUT FILE FORMATS
------------------
  Query file:
    A plain text file containing a single DNA sequence.
    The sequence may span multiple lines. All lines are
    concatenated into one continuous string. Lines
    starting with > are treated as headers and ignored.

  Database file (FASTA format):
    Each sequence begins with a header line starting
    with >. All subsequent lines until the next header
    (or end of file) are concatenated to form the
    sequence. Blank lines are ignored. Both uppercase
    and lowercase letters are accepted and normalized
    to uppercase internally.

ALGORITHMS IMPLEMENTED
----------------------
  1. Longest Common Substring
     - Finds the longest contiguous substring shared
       between two sequences.
     - Similarity score = length of that substring.
     - Time: O(m*n) where m, n are sequence lengths.
     - Space: O(min(m,n)) with rolling rows.

  2. Longest Common Subsequence (LCS)
     - Classic DP algorithm. Finds the longest
       subsequence (not necessarily contiguous) common
       to both sequences.
     - Similarity score = length of the LCS.
     - Time: O(m*n).   Space: O(min(m,n)).

  3. Edit Distance (Levenshtein)
     - Counts the minimum number of insertions,
       deletions, and substitutions to transform one
       sequence into the other.
     - Similarity score = negative edit distance
       (lower distance = more similar).
     - Time: O(m*n).   Space: O(min(m,n)).

  4. Needleman-Wunsch Global Alignment
     - The gold standard for global pairwise alignment.
       Uses a scoring matrix: +1 match, -1 mismatch,
       -2 gap penalty.
     - Similarity score = alignment score.
     - Time: O(m*n).   Space: O(m*n).
