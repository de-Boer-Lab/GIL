import csv

def enumerate_all_sequences(alphabet, seqLen):
    ''' 
    Takes an alphabet and a length of sequences and then returns a list of all
    possible sequences of length seqLen using all characters in alphabet.

        Parameters:
            alphabet (list of str): A list of characters to include in
                                    the generated sequences, e.g. ['A', 'C', 'T', 'G'].
            seqLen (int): Length of the generated sequences.

        Returns:
            mySeqs (list of str): All possible sequences of length seqLen
                                  that can be made from alphabet.
    '''
    if seqLen==1:
        return alphabet
    seqLen=seqLen-1
    parentSeqs = enumerate_all_sequences(alphabet, seqLen)
    mySeqs=[]
    for s in parentSeqs:
        for b in alphabet:
            mySeqs.append(s+b)
    return mySeqs


def reverse_complement(seq):
    '''
    Returns the reverse complement of a DNA sequence.
    
        Parameters:
            seq (str): A DNA sequence.
        
        Returns:
            result (str): The reverse complement of seq.
    '''
    complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    result = []
    for i in range(len(seq)-1, -1, -1):
            result.append(complement_dict[seq[i]])
    return "".join(result)


def hammingdistance(ref, query):
    '''
    Calculates the hamming distance between two sequences. The sequences must be the same length.
    
        Parameters:
            ref (str): First string to compare.
            query (str): Second string to compare.

        Returns:
            count (int): Hamming distance between ref and query.
    '''
    if len(ref) != len(query):
        raise ValueError("ref and query must be same length")

    count = 0

    ref = ref.upper()
    query = query.upper()

    for i in range(0, len(ref), 1):
        if ref[i] == query[i]:
            continue
        else:
            count = count + 1

    return count


def levenshtein(s1, s2):
    '''
    Calculates Levenshtein distance between two sequences. Much slower than the Levenshtein module function.
    Source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python.

        Parameters:
            s1 (str): First sequence for distance calculation.
            s2 (str): Second sequence for distance calculation.

        Returns:
            (int): Levenshtein distance between s1 and s2.
    '''
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]


def read_tsv(filename):
  # Check that extension is .tsv
  if filename[-4:] != ".tsv":
    raise ValueError("Error: file extension must be .tsv")
  with open(filename, "r") as f:
    rd = csv.reader(f, delimiter="\t")
    return [line for line in rd]


def save_seqs_to_tsv(seqs, cols, filename):
    '''
    Saves a list of sequences as a .tsv file with the sequences ordered into
    a matrix with a specified number of columns. 

        Parameters:
            index (str): Sequence of primer index.
            read (str): Valid values are "i7" and "i5". Determines whether the
                        index is for read 1 (i7) or read 2 (i5) primers.
            library_type (str): Valid values are "TruSeq" or "Nextera". Default "TruSeq".
                                Determines whether primer sequences will be for TruSeq or
                                Nextera chemistry.
        Returns:
            primer (str): The indexing primer with index inserted.
    '''
    if len(seqs) % cols != 0:
        raise ValueError("Number of seqs should be a multiple of cols")
    tab_sep_seqs = ["\t".join(seqs[i:i + cols]) for i in range(0, len(seqs), cols)]
    tab_sep_seqs_newline = "\n".join(tab_sep_seqs)
    with open(filename, "w") as f:
        f.write(tab_sep_seqs_newline)