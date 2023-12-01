import sys

def get_fasta_alphabet(file_path):
    """
    Get the alphabet (unique characters) from a multi-fasta file.
    
    Args:
    file_path (str): Path to the fasta file.

    Returns:
    str: A string of sorted unique characters found in the file.
    """
    if not isinstance(file_path, str) or not file_path:
        raise ValueError("Invalid file path provided.")

    alphabet = set()
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('>'):  # Ignore headers
                    alphabet.update(line.strip())
    except Exception as e:
        raise IOError(f"Error reading file: {e}")

    return ''.join(sorted(alphabet))

def verify_sequence_content(file_path):
    """
    Verify the content of the sequences in the input FASTA file.
    This function checks for any characters other than 'A', 'C', 'G', 'T', or 'N'
    in the sequences (not in the headers).
    
    Args:
    file_path (str): Path to the fasta file.

    Returns:
    set: A set of unexpected characters found in the file.
    """
    if not isinstance(file_path, str) or not file_path:
        raise ValueError("Invalid file path provided.")

    unexpected_chars = set()
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('>'):  # Ignore headers
                    for char in line.strip():
                        if char not in {'A', 'C', 'G', 'T', 'N'}:
                            unexpected_chars.add(char)
    except Exception as e:
        raise IOError(f"Error reading file: {e}")

    return unexpected_chars

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    unexpected_chars = verify_sequence_content(file_path)
    print(f"Unexpected Characters in {file_path}: {unexpected_chars}")
    alphabet = get_fasta_alphabet(file_path)
    print(f"Alphabet in {file_path}: {alphabet}")

if __name__ == "__main__":
    main()
