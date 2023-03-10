import random
import copy
import os
from GIL.tools import *


def filter_startingG(seqs, n):
    '''
    Given a list of sequences, returns all sequences that do not start with G.

        Parameters:
            seqs (list of str): A list of DNA sequences.

        Returns:
            passed (list of str): All sequences in seqs that don't start with G.
    '''
    passed = []
    for seq in seqs:
        if seq[:n] != "G" * n:
            passed.append(seq)
    return passed


def filter_endingC(seqs, n):
    '''
    Given a list of sequences, returns all sequences that do not end with C.

        Parameters:
            seqs (list of str): A list of DNA sequences.

        Returns:
            passed (list of str): All sequences in seqs that don't start with C.
    '''
    passed = []
    for seq in seqs:
        if seq[-n:] != "C" * n:
            passed.append(seq)
    return passed


def filter_GC(seqs, minGC=25, maxGC=75):
    '''
    Given a list of sequences, returns all sequences with GC within a desired range.

        Parameters:
            seqs (list of str): A list of DNA sequences.
            minGC (int): GC content must be strictly higher than this value for a sequence to be retained.
            maxGC (int): GC content must be strictly lower than this value for a sequence to be retained.

        Returns:
            passed (list of str): All sequences in seqs that are within the specified GC content range.
    '''
    passed = []
    for seq in seqs:
        G = seq.count("G")
        C = seq.count("C")
        GC_content = ((G+C)/len(seq))*100
        if GC_content > minGC and GC_content < maxGC:
            passed.append(seq)
    return passed


def filter_homopolymer(seqs, max_repeat):
    '''
    Given a list of sequences, returns all sequences without homopolymer repeats greater than a specified length.

        Parameters:
            seqs (list of str): A list of DNA sequences.
            max_repeat (int): Maximum length of homopolymer repeat to retain.

        Returns:
            passed (list of str): All sequences in seqs that contain homopolymer repeats of length
                                  no larger than max_repeat.
    '''
    if max_repeat < 1:
        raise ValueError("max_repeat must be >= 1")
    passed = []
    for seq in seqs:
        found_homopolymer = False
        for i in range(len(seq) - max_repeat):
            if len(set(seq[i:i + max_repeat + 1])) == 1:
                found_homopolymer = True
        if not found_homopolymer:
            passed.append(seq)
    return passed


def test_dinucleotide(seq, max_dinu, cursor):
    '''
    Given a sequence, returns true if the sequence contains a dinucleotide repeat with more than
    a given number of repeats, with dinucleotides defined as starting at the cursor (position 0 or 1).

        Parameters:
            seq (str): A DNA sequence.
            max_dinu (int): Maximum number of dinucleotide repeats to retain.
            cursor (int): Position that dinucleotides are defined from (0 or 1).

        Returns:
            (bool): True if sequence contains > max_dinu dinucleotide repeats.
    '''
    curr_count = 0

    while cursor < len(seq):
        dinu = seq[cursor:cursor + 2]
        if cursor > 1 and dinu == seq[cursor - 2:cursor]:
            curr_count += 1
            if curr_count > max_dinu:
                return False
        else:
            curr_count = 1
        cursor += 2
    return True


def contains_dinucleotide_repeats(seq, max_dinu):
    '''
    Given a sequence, returns true if the sequence contains a dinucleotide repeat with more than
    a given number of repeats.

        Parameters:
            seq (str): A DNA sequence.
            max_dinu (int): Maximum number of dinucleotide repeats to retain.

        Returns:
            (bool): True if sequence contains > max_dinu dinucleotide repeats.
    '''
    if len(seq) < 2*max_dinu:
        print('%s invalid index seq or dinucleotide threshold' % seq)
        return False

    dinu_start_zero = test_dinucleotide(seq, max_dinu, cursor=0)
    dinu_start_one = test_dinucleotide(seq, max_dinu, cursor=1)

    return (dinu_start_zero and dinu_start_one)


def filter_dinucleotide_repeats(seqs, max_dinu):
    '''
    Given a list of sequences, returns all sequences without dinucleotide repeats greater than a specified length.

        Parameters:
            seqs (list of str): A list of DNA sequences.
            max_dinu (int): Maximum number of dinucleotide repeats to retain.

        Returns:
            passed (list of str): All sequences in seqs that do not contain dinucleotide repeats with 
                                  number of repeats larger than max_repeat.
    '''
    passed = []
    for seq in seqs:
        if contains_dinucleotide_repeats(seq, max_dinu):
            passed.append(seq)
    return passed


def self_priming_distance(index, primer_5prime, primer_3prime):
    '''
    Returns the minimum Hamming distance between the reverse complement of 8 bps at the 3' end and of a
    primer and any 8 bp window within the sequence that includes the index sequence. Low hamming distance
    indicates a potential for self-priming.
        
        Parameters:
            index (str): Primer index. This is the generated index which corresponds to the sequence
                         that will be read by the sequencer, NOT the actual index sequence within the primer.
            primer_5prime (str): The 5' end of the primer.
            primer_3prime (str): Sequence of 3' adapter region.
        
        Returns:
            (int): Minimum Hamming distance between the reverse complement of the final 8 bp of the primer
            and any 8 bp window in the primer that includes the index sequence.
    '''

    # Only look at the region +/- 7 bp from the index; we can't change anything about
    # self-priming outside of this region
    index_region = primer_5prime[-7:] + reverse_complement(index) + primer_3prime[0:7]
    distances = []

    # Reverse complement the index because the index sequences that are being filtered correspond to the
    # sequences that will be read by the sequencer, NOT the sequences present in the primer.
    primer_3prime_end_rc = reverse_complement(primer_3prime[-8:])
    for i in range(len(index_region) - 7):
        distances.append(hammingdistance(index_region[i:i + 8], primer_3prime_end_rc))
        
    return min(distances)


def filter_self_priming(indexes, min_dist, primer_5prime, primer_3prime):
    '''
    Given a set of indexes and the primer 5' and 3' sequences,
    returns all indexes that have a low self-priming potential in the index region.
        
        Parameters:
            indexes (list of str): Primer indexes. These are the generated indexes which correspond to the sequences
                                   that will be read by the sequencer, NOT the actual index sequences within the primer.
            min_dist (int): Minimum distance between 3' end of primer and index region of primer (defined by self_priming_distance).
            primer_5prime (str): the sequence of the indexing primer on the 5' side of the index
            primer_3prime (str): the sequence of the indexing primer on the 3' side of the index
        
        Returns:
            passed (list of str): All indexes that do not have a high self-priming potential.
    '''
    if len(primer_5prime) < 8 or len(primer_3prime) < 8:
        raise ValueError("The primer is too short to filter self-priming sequences. Use the --no_filter_self_priming argument "
                         "to skip this filtering step")

    passed = []
    for index in indexes:
        if self_priming_distance(index, primer_5prime, primer_3prime) >= min_dist:
            passed.append(index)
    return passed


def filter_blocklist(seqs, blocklist, dist, dist_fun):
    '''
    Filters out sequences that are too similar to any sequence in a blocklist of sequences.
    
        Parameters:
            seqs (list of str): List of sequences being filtered.
            blocklist (list of str): List of sequences that output sequences must be dissimilar to.
            dist (int): Minimum distance between output seqs and any sequence in blocklist.
            dist_fun: Function used to calculate distance.

        Returns:
            passed (list of str): All sequences in seqs that are at least dist away from any sequence
                                  in blocklist as measured by dist_fun.
    '''
    blocklist = [seq.upper() for seq in blocklist]

    chars = set("".join(blocklist))
    if chars != {'A', 'C', 'G', 'T'}:
        raise ValueError(f"Blocklist indexes must only contain A, C, G, T. These characters are present in your blocklist: {chars}.")

    seq_len = len(seqs[0])
    for seq in blocklist:
        if len(seq) != seq_len:
            raise ValueError(f"Blocklist indexes must be the same length as generated sequences.")

    passed = []
    for seq in seqs:
        dissimilar = True
        for index in blocklist:
            if dist_fun(seq, index) < dist:
                dissimilar = False
                break
        if dissimilar:
            passed.append(seq)
    return passed


def select_dissimilar_sequences(seqs, dist, dist_fun):
    '''
    Selects a set of sequences randomly from seqs such that all pairwise distances between
    sequences in the set meet a minimum threshold. The size of the output set will vary depending
    on the random sampling of seqs.
    
        Parameters:
            seqs (list of str):  Set of sequences to sample from.
            max_seqs (int): Maximum number of sequences to return (stop searching after this number is reached).
            dist (int): Minimum pairwise distance between all sequences in the output set
            dist_fun (function): A function that takes as parameters two sequences and returns
                                 a distance between the sequences (default Levenshtein distance).
        
        Returns:
            final (list of str): A list of sequences sampled from seqs such that all pairwise distances between
            between sequences in the set meet a minimum threshold (dist).
            
            '''
    seqs_remaining = copy.deepcopy(seqs)
    final = []
    while len(seqs_remaining) > 0:
        seq = random.sample(seqs_remaining, 1)[0]
        final.append(seq)
        seqs_remaining.remove(seq)
        i=0
        while i < len(seqs_remaining):
            if dist_fun(seqs_remaining[i], seq) < dist:
                del seqs_remaining[i]
            else: # only increment if we didn't just delete one
                i += 1
    return final
    

def colour_balanced(seqs, channel_2 = True, channel_4 = True):
    '''
    Determines if sequences meet colour diversity requirements for 2 and/or 4 channel sequencers.
    
        Parameters:
            seqs (list of str): List of DNA sequences to check colour compatibility.
                                All sequences must be same length.
            channel_2 (bool): Set to true to check diversity requirements for 2 channel sequencer.
            channel_4 (bool): Set to true to check diversity requirements for 4 channel sequencer.
        
        Returns:
            (bool): True if seqs are colour balanced, False if not colour balanced.
    
    '''
    if channel_2:
        for i in range(len(seqs[0])): # For each postion in sequence
            saw_C = False #Set C found to False
            saw_T = False #Set T found to False
            for j, seq in enumerate(seqs): # look at  position i in each sequence
                if seq[i] == "A": #If A/T/C is at postion i move to next position
                    break
                elif seq[i] == "C":
                    saw_C = True
                elif seq[i] == "T":
                    saw_T = True
                if saw_C and saw_T:
                    break
                if j == len(seqs) - 1: 
                    return False

    if channel_4:
        for i in range(len(seqs[0])): # For each position in sequence
            a_C = False #Set A or C found to False
            g_T = False #Set G or T found to False
            for seq in seqs: # look at  position i in each sequence
                if a_C == True and g_T == True: #If we've already found both a 'A' or 'C' and a 'G' or 'T' go to next position
                    break
                if seq[i] in ['A','C'] and a_C == False: #If position i in sequence [seq] is an 'A' or 'C' set a_C to True (found)
                    a_C = True
                if seq[i] in ['G','T'] and g_T == False: #If position i in sequence [seq] is an 'G' or 'T' set g_T to True (found)
                    g_T = True
            if not a_C or not g_T: # If we haven't found an A/C or G/T at position i, return false
                return False

    return True


def generate_colour_balanced_indices(seqs, group_size, max_iterations, channel_2 = True, channel_4 = True):
    '''
    Returns a list of lists of indexes that are colour balanced. Each inner list of length group_size
    is colour balanced.

        Parameters:
            seqs (list of str): List of DNA sequences to check colour compatibility.
                                All sequences must be same length.
            max_iterations (int): Maximum number of times a random sample of sequences will
                                  be taken to check for colour balance before giving up
                                  and returning the generated sequences.
            channel_2 (bool): Set to True to check colour balance for 2 channel sequencer.
            channel_4 (bool): Set to True to check colour balance for 4 channel sequencer.

        Returns:
            (list of lists of str): A list of lists of sequences where the sequences in each inner
                                    list of length group_size are all colour balanced.
    '''
    counter = 0
    colour_balanced_groups = []
    seqs = copy.deepcopy(seqs)

    while counter < max_iterations:
        if len(seqs) < group_size:
            # print("Ran out of sequences")
            break

        sampled_seqs = random.sample(seqs, group_size)

        if colour_balanced(sampled_seqs):
            colour_balanced_groups.append(sampled_seqs)
            for seq in sampled_seqs:
                seqs.remove(seq)
            counter = 0
        else:
            counter += 1

    return colour_balanced_groups


def create_primer_plates(i7_lists, i5_lists, i5_start, i5_end, i7_start, i7_end, dir, plate_name):
    '''
    Makes 96-well plates from a list of lists of i7 and i5 indexes. The length of the inner lists represents the
    smallest number of indexes that can be chosen while maintaining colour balance (default 4). 
    '''

    i7s = [index for sublist in i7_lists for index in sublist]
    i5s = [index for sublist in i5_lists for index in sublist]

    rc_i7s = [reverse_complement(index) for index in i7s]
    rc_i5s = [reverse_complement(index) for index in i5s]

    i7_primers = [i7_start + index + i7_end for index in rc_i7s]
    i5_primers = [i5_start + index + i5_end for index in rc_i5s]

    num_plates = min(len(i7s) // 96, len(i5s) // 96)

    if num_plates < 1:
        raise ValueError("Not enough indexes to create a 96 well plate")

    os.makedirs(f"{dir}/Plates/Indexes", exist_ok=True)
    os.makedirs(f"{dir}/Plates/Indexes/i7", exist_ok=True)
    os.makedirs(f"{dir}/Plates/Indexes/i5", exist_ok=True)
    os.makedirs(f"{dir}/Plates/Primers", exist_ok=True)

    for i in range(num_plates):
        save_seqs_to_tsv(rc_i7s[i*96:(i + 1)*96], 12, f"{dir}/Plates/Indexes/i7/{plate_name}_i7_Indexes_Plate_{i + 1}.tsv")
        save_seqs_to_tsv(i7_primers[i*96:(i + 1)*96], 12, f"{dir}/Plates/Primers/{plate_name}_i7_Primers_Plate_{i + 1}.tsv")
        save_seqs_to_tsv(rc_i5s[i*96:(i + 1)*96], 12, f"{dir}/Plates/Indexes/i5/{plate_name}_i5_Indexes_Plate_{i + 1}.tsv")
        save_seqs_to_tsv(i5_primers[i*96:(i + 1)*96], 12, f"{dir}/Plates/Primers/{plate_name}_i5_Primers_Plate_{i + 1}.tsv")


def make_order_sheet(primers_file, sample_sheet_name, primer_name, dir, company="IDT", mod=True):

    primers = read_tsv(primers_file)

    if company == "IDT":
        primers_sheet = ["Well Position,Name,Sequence"]
        for i, letter in enumerate(["A", "B", "C", "D", "E", "F", "G", "H"]):
            for j, primer in enumerate(primers[i]):
                if mod:
                    primer = primer[:-1] + "*" + primer[-1]
                primers_sheet.append(f'{letter}{j + 1},{letter}{j + 1}_{primer_name},{primer}')

    elif company == "Thermo":
        primers_sheet = ["Row,Column,Name,Sequence"]
        mod_dict = {'A': 'F', 'C': 'O', 'G': 'E', 'T': 'Z'}
        for i, letter in enumerate(["A", "B", "C", "D", "E", "F", "G", "H"]):
            for j, primer in enumerate(primers[i]):
                if mod:
                    primer = list(primer)
                    primer[-2] = mod_dict[primer[-2]]
                    primer = "".join(primer)
                primers_sheet.append(f'{letter},{j + 1},{letter}{j + 1}_{primer_name},{primer}')

    elif company == "Sigma":
        primers_sheet = ["Row,Column,Name,Sequence"]
        for i in range(12):
            for j, letter in enumerate(["A", "B", "C", "D", "E", "F", "G", "H"]):
                primer = primers[j][i]
                if mod:
                    primer = primer[:-1] + "*" + primer[-1]
                primers_sheet.append(f'{letter},{i + 1},{letter}{i + 1}_{primer_name},{primer}')

    else:
        raise ValueError("Company must be IDT, Thermo, or Sigma")

    os.makedirs(f"{dir}/Order_Sheets", exist_ok=True)

    with open(f"{dir}/Order_Sheets/{sample_sheet_name}_order_sheet.csv", 'w') as f:
      f.write("\n".join(primers_sheet))