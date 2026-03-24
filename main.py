"""
main.py
Entry point for the DNA Sequence Matching project.

Presents a menu-driven interface that:
  1. Asks for the query file and database file paths.
  2. Parses both files.
  3. Lets the user choose a similarity algorithm.
  4. Searches the database for the most similar sequence.
  5. Reports the result and loops back to the menu.
"""

import time
from file_handler import parse_fasta, parse_query
from algorithms import ALGORITHMS


def print_banner():
    """Print the program header."""
    print()
    print("=" * 50)
    print("  DNA Sequence Matching — CSCI 311 Spring 2026")
    print("=" * 50)
    print()


def print_menu():
    """Print the algorithm selection menu."""
    print()
    print("=" * 50)
    print("  Select a Similarity Algorithm")
    print("=" * 50)
    for key, (name, _) in ALGORITHMS.items():
        print(f"   {key}. {name}")
    print(f"   {len(ALGORITHMS) + 1}. Exit")
    print("=" * 50)


def get_file_paths() -> tuple[str, str]:
    """
    Prompt the user for the query file and database file paths.
    Keeps asking until valid files are provided.
    """
    while True:
        query_path = input("Enter path to query sequence file: ").strip()
        try:
            query_seq = parse_query(query_path)
            print(f"  -> Query loaded: {len(query_seq)} nucleotides.")
            break
        except FileNotFoundError:
            print(f"  [Error] File not found: '{query_path}'. Try again.")
        except ValueError as e:
            print(f"  [Error] {e}. Try again.")

    while True:
        db_path = input("Enter path to database sequence file: ").strip()
        try:
            db_seqs = parse_fasta(db_path)
            print(f"  -> Database loaded: {len(db_seqs)} sequence(s).")
            break
        except FileNotFoundError:
            print(f"  [Error] File not found: '{db_path}'. Try again.")
        except ValueError as e:
            print(f"  [Error] {e}. Try again.")

    return query_seq, db_seqs


def find_most_similar(query, db_seqs, algo_func, algo_name):
    """
    Search the database for the sequence most similar to the query.

    Implements the FindMostSimilarSeq pseudocode from the assignment:
      For each sequence in D, compute similarity with the query.
      Track and return the best match.

    Parameters
    ----------
    query : str
        The query DNA sequence.
    db_seqs : list of (str, str)
        The database: list of (header, sequence) tuples.
    algo_func : callable
        The similarity function to use.
    algo_name : str
        Human-readable name of the algorithm (for display).

    Returns
    -------
    best_header : str
        Header of the most similar sequence.
    best_score : int/float
        The similarity score of the best match.
    elapsed : float
        Time taken in seconds.
    """
    print(f"\n  Searching with: {algo_name}")
    print(f"  Comparing against {len(db_seqs)} sequence(s)...\n")

    best_score = float("-inf")
    best_header = None
    start_time = time.time()

    for i, (header, seq) in enumerate(db_seqs, start=1):
        # Extract a short name from the header for progress display
        short_name = header[1:].split()[0] if len(header) > 1 else "???"
        print(f"    [{i}/{len(db_seqs)}] {short_name} ", end="", flush=True)

        score = algo_func(seq, query)

        print(f"... score = {score}")

        if score > best_score:
            best_score = score
            best_header = header

    elapsed = time.time() - start_time

    return best_header, best_score, elapsed


def main():
    """Main program loop."""
    print_banner()

    # Step 1: Get file paths and parse files
    query_seq, db_seqs = get_file_paths()

    # Step 2: Algorithm selection loop
    exit_choice = len(ALGORITHMS) + 1

    while True:
        print_menu()

        try:
            choice = int(input(f"  Select an algorithm (1-{exit_choice}): ").strip())
        except ValueError:
            print("  [Error] Please enter a number.")
            continue

        if choice == exit_choice:
            print("\n  Goodbye!\n")
            break

        if choice not in ALGORITHMS:
            print(f"  [Error] Invalid choice. Enter a number between 1 and {exit_choice}.")
            continue

        algo_name, algo_func = ALGORITHMS[choice]

        # Step 3: Run the search
        best_header, best_score, elapsed = find_most_similar(
            query_seq, db_seqs, algo_func, algo_name
        )

        # Step 4: Report results
        print()
        print("=" * 50)
        print("  RESULT")
        print("=" * 50)
        print(f"  Algorithm : {algo_name}")
        print(f"  Most similar sequence:")
        print(f"    {best_header}")
        print(f"  Score     : {best_score}")
        print(f"  Time      : {elapsed:.2f} seconds")
        print("=" * 50)


if __name__ == "__main__":
    main()
