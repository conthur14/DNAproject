"""
file_handler.py
Handles reading and parsing of DNA sequence files.

Supports:
  - FASTA format (database files with >header lines)
  - Plain text query files (raw sequence, optionally with a header)
"""


def parse_fasta(filepath: str) -> list[tuple[str, str]]:
    """
    Parse a FASTA-formatted file into a list of (header, sequence) tuples.

    FASTA format rules:
      - Header lines start with '>'
      - Sequence lines follow the header and are concatenated
      - Blank lines are ignored
      - Letters are normalized to uppercase

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.

    Returns
    -------
    list of (str, str)
        Each element is (header_string, sequence_string).
        The header includes the '>' character.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file contains no valid sequences.
    """
    sequences = []
    current_header = None
    current_seq_parts = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Skip blank lines
            if not line:
                continue

            # Header line
            if line.startswith(">"):
                # Save the previous sequence if there is one
                if current_header is not None:
                    full_seq = "".join(current_seq_parts).upper()
                    sequences.append((current_header, full_seq))

                current_header = line
                current_seq_parts = []
            else:
                # Sequence data line — accumulate it
                current_seq_parts.append(line)

    # Don't forget the last sequence in the file
    if current_header is not None:
        full_seq = "".join(current_seq_parts).upper()
        sequences.append((current_header, full_seq))

    if not sequences:
        raise ValueError(f"No valid FASTA sequences found in '{filepath}'.")

    return sequences


def parse_query(filepath: str) -> str:
    """
    Parse a query file into a single DNA sequence string.

    The query file may be:
      - A plain text file with the sequence on one or more lines
      - A FASTA file with a single header + sequence

    In either case, header lines (starting with '>') are skipped,
    and all remaining lines are concatenated and uppercased.

    Parameters
    ----------
    filepath : str
        Path to the query file.

    Returns
    -------
    str
        The query sequence as a single uppercase string.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If no sequence data is found.
    """
    seq_parts = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue  # skip header lines
            seq_parts.append(line)

    query = "".join(seq_parts).upper()

    if not query:
        raise ValueError(f"No sequence data found in '{filepath}'.")

    return query


# ---------------------------------------------------------------------------
# Quick self-test: run this file directly to verify parsing works
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python3 file_handler.py <fasta_file> [query_file]")
        print("  Tests the parser on the given file(s).")
        sys.exit(1)

    # Test FASTA parsing
    fasta_path = sys.argv[1]
    print(f"Parsing FASTA file: {fasta_path}")
    seqs = parse_fasta(fasta_path)
    print(f"  Found {len(seqs)} sequence(s):\n")
    for header, seq in seqs:
        print(f"  {header}")
        print(f"    Length: {len(seq)} nucleotides")
        print(f"    Preview: {seq[:60]}...\n")

    # Test query parsing if a second file is provided
    if len(sys.argv) >= 3:
        query_path = sys.argv[2]
        print(f"Parsing query file: {query_path}")
        q = parse_query(query_path)
        print(f"  Length: {len(q)} nucleotides")
        print(f"  Preview: {q[:60]}...")
