"""
algorithms.py
Implements sequence similarity algorithms for DNA comparison.

Each function takes two sequences (s, t) and returns a numeric
similarity score. Higher scores = more similar.

Algorithms:
  1. Longest Common Substring
  2. Longest Common Subsequence (LCS)
  3. Edit Distance (negated, so higher = more similar)
  4. Needleman-Wunsch Global Alignment
"""


def longest_common_substring(s: str, t: str) -> int:
    """
    Find the length of the longest contiguous substring common to s and t.

    Uses a 2-row rolling DP table to save space.

    Time:  O(m * n)
    Space: O(min(m, n))

    Parameters
    ----------
    s, t : str
        The two DNA sequences to compare.

    Returns
    -------
    int
        Length of the longest common substring (higher = more similar).
    """
    # Make sure the shorter string is used for the column dimension
    if len(s) < len(t):
        s, t = t, s

    m, n = len(s), len(t)
    # prev and curr are the two rows of the DP table
    prev = [0] * (n + 1)
    curr = [0] * (n + 1)
    best = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                curr[j] = prev[j - 1] + 1
                if curr[j] > best:
                    best = curr[j]
            else:
                curr[j] = 0
        # Swap rows: curr becomes prev, and we reuse the old prev
        prev, curr = curr, prev
        # Reset curr for next iteration
        for j in range(n + 1):
            curr[j] = 0

    return best


def longest_common_subsequence(s: str, t: str) -> int:
    """
    Compute the length of the Longest Common Subsequence of s and t.

    Unlike substring, the matching characters do not need to be
    contiguous — they just need to appear in the same relative order.

    Standard bottom-up DP approach.

    Time:  O(m * n)
    Space: O(m * n)  [could be reduced to O(min(m,n)) if only length needed]

    Parameters
    ----------
    s, t : str
        The two DNA sequences to compare.

    Returns
    -------
    int
        Length of the LCS (higher = more similar).
    """
    m, n = len(s), len(t)

    # dp[i][j] = LCS length of s[0..i-1] and t[0..j-1]
    # Use a 2-row optimization
    prev = [0] * (n + 1)
    curr = [0] * (n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                curr[j] = prev[j - 1] + 1
            else:
                curr[j] = max(prev[j], curr[j - 1])
        prev, curr = curr, [0] * (n + 1)

    return prev[n]


def edit_distance(s: str, t: str) -> int:
    """
    Compute the edit distance (Levenshtein distance) between s and t,
    then return its negation so that MORE similar sequences score HIGHER.

    Edit operations: insertion, deletion, substitution (each costs 1).

    Time:  O(m * n)
    Space: O(min(m, n))

    Parameters
    ----------
    s, t : str
        The two DNA sequences to compare.

    Returns
    -------
    int
        Negative edit distance (closer to 0 = more similar).
    """
    # Ensure the shorter string is used for the inner loop
    if len(s) < len(t):
        s, t = t, s

    m, n = len(s), len(t)

    # Only need two rows
    prev = list(range(n + 1))  # base case: transforming "" into t[0..j-1]
    curr = [0] * (n + 1)

    for i in range(1, m + 1):
        curr[0] = i  # transforming s[0..i-1] into ""
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                curr[j] = prev[j - 1]  # no operation needed
            else:
                curr[j] = 1 + min(
                    prev[j],      # deletion
                    curr[j - 1],  # insertion
                    prev[j - 1]   # substitution
                )
        prev, curr = curr, [0] * (n + 1)

    # Negate so that higher = more similar
    return -prev[n]


def needleman_wunsch(s: str, t: str) -> int:
    """
    Needleman-Wunsch global alignment score.

    Scoring scheme:
        Match:    +1
        Mismatch: -1
        Gap:      -2

    This is the standard bioinformatics global alignment algorithm.

    Time:  O(m * n)
    Space: O(m * n)  [full table needed for potential traceback]

    Parameters
    ----------
    s, t : str
        The two DNA sequences to compare.

    Returns
    -------
    int
        The alignment score (higher = more similar).
    """
    MATCH = 1
    MISMATCH = -1
    GAP = -2

    m, n = len(s), len(t)

    # Build the full DP table
    # dp[i][j] = best alignment score for s[0..i-1] vs t[0..j-1]
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Base cases: aligning against an empty sequence = all gaps
    for i in range(1, m + 1):
        dp[i][0] = i * GAP
    for j in range(1, n + 1):
        dp[0][j] = j * GAP

    # Fill the table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                diag = dp[i - 1][j - 1] + MATCH
            else:
                diag = dp[i - 1][j - 1] + MISMATCH

            up   = dp[i - 1][j] + GAP      # gap in t
            left = dp[i][j - 1] + GAP      # gap in s

            dp[i][j] = max(diag, up, left)

    return dp[m][n]


# ---------------------------------------------------------------------------
# Registry: maps algorithm names to their functions and description
# ---------------------------------------------------------------------------
ALGORITHMS = {
    1: ("Longest Common Substring",           longest_common_substring),
    2: ("Longest Common Subsequence (LCS)",   longest_common_subsequence),
    3: ("Edit Distance",                      edit_distance),
    4: ("Needleman-Wunsch Alignment",         needleman_wunsch),
}
