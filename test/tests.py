import os
import pathlib
import random
import re
import unittest
from GIL.pipeline_components import *
from GIL.tools import *

levenshtein_imported = False
try:
 import Levenshtein
 levenshtein_dist = Levenshtein.distance
 levenshtein_imported = True
except ModuleNotFoundError:
    print("Levenshtein module (https://github.com/maxbachmann/Levenshtein) not found. "
            "A much slower implementation of Levenshtein distance will be used.")
    levenshtein_dist = levenshtein

class TestFilters(unittest.TestCase):
    def test_filter_startingG(self):
        sequences = ["GTATAATATAT","GCGCGCGCGCGC"]
        sequences2 = ["GGATAATATCC","GCGCGCGCGCGC"]
        self.assertEqual(len(filter_startingG(sequences, 1)), 0, msg="FAILED G START")
        self.assertEqual(len(filter_startingG(sequences2, 2)), 1, msg="FAILED GG START")
        self.assertEqual(len(filter_startingG(["GREAT","GREY","GOOSE"], 1)), 0, msg="filter_startingG did not eliminate all G-starting words")

    def test_filter_endingC(self):
        endingC_test_seqs = ['ATACGACC','ATACGACG']
        sequences = ["GGATAATATCC","GCGCGCGCGCGC"]
        self.assertEqual(len(filter_endingC(sequences, 1)), 0, msg="FAILED C END")
        self.assertEqual(len(filter_endingC(sequences, 2)), 1, msg="FAILED CC END")
        self.assertEqual(filter_endingC(endingC_test_seqs, 1), ['ATACGACG'])
        
    def test_filter_GC(self):
        sequences = ["GTATAATATAT","GCGCGCGCGCGC"]
        self.assertEqual(len(filter_GC(sequences)), 0, msg="FAILED GC CONTENT")
        self.assertEqual(len(filter_GC(["GGGGA","CCCCA","TTTTG","GATCG","GGTCC"])), 1, msg="GC did not eliminate all high GC sequences")
        
        # seqs = filter_GC(seqs)
        # for seq in seqs:
        #     G = seq.count("G")
        #     C = seq.count("C")
        #     GC_content = ((G+C)/len(seq))*100

        # assert GC_content > 25 or GC_content < 75, "GC content not following the range"
    
        self.assertEqual(len(filter_GC(['ATGCGCGC', 'ATTTACGT','GATGCGCGC','ATTTACGTG'])), 1) # Boundary testing. Just at or exceed 75% or just at or below 25%
        self.assertEqual(len(filter_GC(['AGGGCGCGC', 'ATTTACGTA'])), 0) # Boundary testing. Just below 75% or just above 25%

    def test_filter_homopolymer(self):
        self.assertEqual(len(filter_homopolymer(["AAATTTTTAAA"], 4)), 0, msg="contains_homopolymer did not find the T5 in AAATTTTTAAA")
        self.assertEqual(len(filter_homopolymer(['AAATAACA'], 3)), 1)
        self.assertEqual(len(filter_homopolymer(['AAATAACA'], 2)), 0)
        self.assertEqual(len(filter_homopolymer(['AAATAACAAAA'], 19)), 1)

    def test_dinucleotide_repeats(self):
        self.assertTrue(contains_dinucleotide_repeats("ATGCATAT", 2), msg="indexingfilter_dinucleotide_repeats incorrectly predicted repeat greater than 2")
        self.assertFalse(contains_dinucleotide_repeats("ATATAT", 2), msg="indexingfilter_dinucleotide_repeats didn't find repeat greater than 2")
        self.assertEqual(len(filter_dinucleotide_repeats(["ATCGCGCG","ACGCGCGA","GCGCGCTA","AAAAAATG"], 2)), 0, msg="dinucleotide_repeats_filter did not filter out all >2xdi repeats")
        
        self.assertTrue(contains_dinucleotide_repeats('ACTGACTG', 2)) # basic
        self.assertFalse(contains_dinucleotide_repeats('ACACACTG', 2)) # basic, repeat starting at the first pos
        self.assertTrue(contains_dinucleotide_repeats('ATGATATC', 2)) # not tandom repeat
        self.assertFalse(contains_dinucleotide_repeats('', 3)) # bad input
        self.assertFalse(contains_dinucleotide_repeats('ACTGTGTGAGTG', 2)) # tandom repeat in middle
        self.assertTrue(contains_dinucleotide_repeats('AAAAGGTG', 2)) #homopolymer
        self.assertFalse(contains_dinucleotide_repeats('AATATATGTG', 2)) # repeat starting at the second pos
        self.assertFalse(contains_dinucleotide_repeats('ATGCGCGC', 2)) #repeat ends at the last pos

        return True

    def test_filter_blocklist(self):
        seqs = ['AAAAAAAA', 'CCCCCCCC', 'GGGGGGGG', 'TTTTTTTT', 'ATATATAT']
        blocklist = ['AAAAAAAA', 'CCCCGGGG', 'TTATTGTA', 'GGAGGGTG', 'TATATATA']
        blocklist_lower = ['aaaaaaaa', 'ccccgggg', 'ttattgta', 'ggagggtg', 'tatatata']
        filtered = filter_blocklist(seqs, blocklist, dist=3, dist_fun=levenshtein_dist)
        self.assertEqual(filtered, ['CCCCCCCC', 'TTTTTTTT'])
        filtered = filter_blocklist(seqs, blocklist_lower, dist=3, dist_fun=levenshtein_dist)
        self.assertEqual(filtered, ['CCCCCCCC', 'TTTTTTTT'])

        # check testing for non-DNA characters in blocklist
        blocklist_n = ['aaaaaaan', 'ccccgggg', 'ttattgta', 'ggagggtg', 'tatatata']
        with self.assertRaises(ValueError):
            filter_blocklist(seqs, blocklist_n, dist=3, dist_fun=levenshtein_dist)

        # check testing for index length in blocklist
        blocklist_long = ['aaaaaaaa', 'ccccgggg', 'ttattgta', 'ggagggtg', 'tatatatat']
        with self.assertRaises(ValueError):
            filter_blocklist(seqs, blocklist_long, dist=3, dist_fun=levenshtein_dist)

    def test_filter_self_priming(self):
        pstart_i5 = 'AATGATACGGCGACCACCGAGATCTACAC'
        pend_i5 = 'TCGTCGGCAGCGTC'
        pstart_i7 = 'CAAGCAGAAGACGGCATACGAGAT'
        pend_i7 = 'GTCTCGTGGGCTCGG'
        i7_test_seq_set = ['GGGCTCGG', 'GGACTCGG', 'GGGCTCTA']
        i5_test_seq_set = ['GCAGCGTC', 'GCCGCGTC', 'GCAGCACC']
        
        self.assertEqual(filter_self_priming(i7_test_seq_set, 3, primer_5prime=pstart_i7, primer_3prime=pend_i7), [])
        self.assertEqual(filter_self_priming(i5_test_seq_set, 3, primer_5prime=pstart_i5, primer_3prime=pend_i5), [])

    def test_self_priming_distance(self):
        pstart_i5 = 'AATGATACGGCGACCACCGAGATCTACAC'
        pend_i5 = 'TCGTCGGCAGCGTC'
        pstart_i7 = 'CAAGCAGAAGACGGCATACGAGAT'
        pend_i7 = 'GTCTCGTGGGCTCGG'
        i5_testseq = pend_i5[-8:] #'GACGCTGC'
        i7_testseq = pend_i7[-8:] #'CCGAGCCC'
        
        self.assertEqual(self_priming_distance(i5_testseq, pstart_i5, pend_i5), 0)
        self.assertEqual(self_priming_distance(i7_testseq, pstart_i7, pend_i7), 0)
        self.assertEqual(self_priming_distance('AAAAAAAA', pstart_i7, pend_i7), 4)
        self.assertEqual(self_priming_distance('AAAAAAAG', pstart_i7, pend_i7), 3)

        self.assertEqual(len(filter_self_priming([i5_testseq], 1, primer_5prime=pstart_i5, primer_3prime=pend_i5)), 0)
        self.assertEqual(len(filter_self_priming([i7_testseq], 1, primer_5prime=pstart_i7, primer_3prime=pend_i7)), 0)
        self.assertEqual(len(filter_self_priming(['CCGAGCCC'], 1, primer_5prime=pstart_i7, primer_3prime=pend_i7)), 1)

class TestTools(unittest.TestCase):
    def test_hammingdistance(self):
        self.assertEqual(hammingdistance("ATGCCGATTAC", "ATGCTTTTTAC"), 3, msg="hammdist got the hammdist wrong between ATGCCGATTAC and ATGCTTTTTAC (should be 3)")
        self.assertGreaterEqual(len(select_dissimilar_sequences(["ATGCCGATTAC", "ATGCTTTTTAC","ATGCCGTTTAC"], dist=3, dist_fun=hammingdistance)), 1, msg="hamming_filter found fewer than 1 sequence (should be 2 or 3)")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("TAAAAAAA"), "TTTTTTTA")
        self.assertEqual(reverse_complement("GC"), "GC")
        self.assertEqual(reverse_complement("GAATTC"), "GAATTC")


    def test_enumerateSeqs(self):
        self.assertEqual(len(set(enumerate_all_sequences(["A", "T", "G", "C"], 8))), 4**8, msg="enumerateSeqs produced the wrong number of unique sequences")

    def test_Levenshtein(self):
        if levenshtein_imported:
        # check Levenshtein implementations
            seqs1 = [random.choices(['A', 'T', 'G', 'C'], k=50) for i in range(500)]
            seqs2 = [random.choices(['A', 'T', 'G', 'C'], k=50) for i in range(500)]
            lev_package_distance = [Levenshtein.distance(seq1, seq2) for seq1, seq2 in zip(seqs1, seqs2)]
            lev_py_distance = [levenshtein(seq1, seq2) for seq1, seq2 in zip(seqs1, seqs2)]
            self.assertEqual(sum((i != j) for i , j in zip(lev_package_distance, lev_py_distance)), 0, msg="Difference between Levenshtein package output and python Levenshtein implementation should be 0")
        else:
            print("Can't compare levenshtein implementations because Levenshtein package could not be imported")

class TestComponents(unittest.TestCase):
    def test_select_dissimilar_sequences(self):
        #Hamming distance
        test_set1 = ['CCGAGCCCA', 'CCGAGGCCA', 'CCGTGGCCA', 'CCGTGACCA' ]
        test_set2 = ['ATGCATCG', "GCTAGCTA", "TACGTACG", "CGATCGAT"]
        self.assertEqual(len(select_dissimilar_sequences(test_set1, dist=3, dist_fun= hammingdistance)), 1)
        self.assertEqual(len(select_dissimilar_sequences(test_set2, dist=3, dist_fun= hammingdistance)), 4)
        
        #Levenshtein distance
        self.assertEqual(len(select_dissimilar_sequences(['ATAT', 'TATA','GGGG', 'GGGG'], dist = 5, dist_fun=levenshtein_dist)), 1)
        self.assertEqual(len(select_dissimilar_sequences(['ATAT', 'TATA','GGGG', 'GGGG'], dist = 3, dist_fun=levenshtein_dist)), 2)

        all_seqs = enumerate_all_sequences(["A", "T", "G", "C"], 4)
        self.assertEqual(len(select_dissimilar_sequences(all_seqs, dist = 0, dist_fun=levenshtein_dist)), len(all_seqs))
        self.assertNotEqual((len(select_dissimilar_sequences(all_seqs, dist = 2, dist_fun=levenshtein_dist))), len(all_seqs))

    def test_colour_balanced(self):
        # 2 channels
        seqs_compat_a = ['AAAAAAAA', 'AAAAAAAA','AAAAAAAA', 'AAAAAAAA']
        seqs_compat_ct = ['CCCCCCCC', 'TTTTTTTT','TTTTTTTT', 'CCCCCCCC']
        seqs_compat_mix = ['CCCCCCCC', 'TTTTTTTT','CCCCCAAA', 'TTTTTAAA']
        seqs_compat_first_A = ['ATGCATGC', 'GCATGCAT', 'GTGCATGC', 'GCATGCAT']
        seqs_not_compat_first_C = ['CTGCATGC', 'CCATGCAT', 'CTGCATGC', 'CCATGCAT']
        seqs_not_compat_first_T = ['TTGCATGC', 'TCATGCAT', 'TTGCATGC', 'TCATGCAT']
        seqs_not_compat_first_G = ['GTGCATGC', 'GCATGCAT', 'GTGCATGC', 'GCATGCAT']
        seqs_compat_mid_A = ['TGCAATGA', 'CATGACAT','TGCAATGA', 'CATGACAT']
        seqs_not_compat_mid_C = ['TGCACTGA', 'CATGCCAT','TGCACTGA', 'CATGCCAT']
        seqs_not_compat_mid_T = ['TGCATTGA', 'CATGTCAT','TGCATTGA', 'CATGTCAT']
        seqs_not_compat_mid_G = ['TGCAGTGA', 'CATGGCAT','TGCAGTGA', 'CATGGCAT']
        seqs_compat_last_A = ['TGCATGCA', 'CATGCATA','TGCATGCA', 'CATGCATA']
        seqs_not_compat_last_C = ['TGCATGCC', 'CATGCATC','TGCATGCC', 'CATGCATC']
        seqs_not_compat_last_G = ['TGCATGCG', 'CATGCATG','TGCATGCG', 'CATGCATG']
        seqs_not_compat_last_T = ['TGCATGCT', 'CATGCATT','TGCATGCT', 'CATGCATT']

        self.assertTrue(colour_balanced(seqs_compat_a, channel_2 = True, channel_4 = False))
        self.assertTrue(colour_balanced(seqs_compat_ct, channel_2 = True, channel_4 = False))
        self.assertTrue(colour_balanced(seqs_compat_mix, channel_2 = True, channel_4 = False))
        self.assertTrue(colour_balanced(seqs_compat_first_A, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_first_C, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_first_T, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_first_G, channel_2 = True, channel_4 = False))
        self.assertTrue(colour_balanced(seqs_compat_mid_A, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_C, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_T, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_G, channel_2 = True, channel_4 = False))
        self.assertTrue(colour_balanced(seqs_compat_last_A, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_last_C, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_last_G, channel_2 = True, channel_4 = False))
        self.assertFalse(colour_balanced(seqs_not_compat_last_T, channel_2 = True, channel_4 = False))
        
        # 4 channels

        seqs_compat_at = ['AAAAAAAA', 'TTTTTTTT','AAAAAAAA', 'TTTTTTTT']
        seqs_compat_gc = ['AAAAAAAA', 'TTTTTTTT','CCCCCCCC', 'CCCCCCCC']
        seqs_compat_mix = ['GGGGGGGG', 'TTTTTTTT','AAAAAAAA', 'CCCCCCCC']
        seqs_not_compat_first_A = ['AATCGATG', 'ACGTACGT','AATCGATG', 'ACGTACGT']
        seqs_not_compat_first_T = ['TATCGATG', 'TCGTACGT','TATCGATG', 'TCGTACGT']
        seqs_not_compat_first_AT = ['AATCGATG', 'TCGTACGT','AATCGATG', 'ACGTACGT']
        seqs_not_compat_first_G = ['GATCGATG', 'GCGTACGT','GATCGATG', 'GCGTACGT']
        seqs_not_compat_first_C = ['CATCGATG', 'CCGTACGT','CATCGATG', 'CCGTACGT']
        seqs_not_compat_first_GC = ['CATCGATG', 'CCGTACGT','CATCGATG', 'GCGTACGT']
        seqs_not_compat_mid_A = ['ATGCAATG', 'CGTAACGT','ATGCAATG', 'CGTAACGT']
        seqs_not_compat_mid_T = ['ATGCTATG', 'CGTATCGT','ATGCTATG', 'CGTATCGT']
        seqs_not_compat_mid_AT = ['ATGCAATG', 'CGTATCGT','ATGCAATG', 'CGTAACGT']
        seqs_not_compat_mid_G = ['ATGCGATG', 'CGTAGCGT','ATGCGATG', 'CGTAGCGT']
        seqs_not_compat_mid_C = ['ATGCCATG', 'CGTACCGT','ATGCCATG', 'CGTACCGT']
        seqs_not_compat_mid_GC = ['ATGCGATG', 'CGTACCGT','ATGCCATG', 'CGTAGCGT']
        seqs_not_compat_last_A = ['ATGCATGA', 'CGTACGTA','ATGCATGA', 'CGTACGTA']
        seqs_not_compat_last_T = ['ATGCATGT', 'CGTACGTT','ATGCATGT', 'CGTACGTT']
        seqs_not_compat_last_AT = ['ATGCATGA', 'CGTACGTT','ATGCATGA', 'CGTACGTA']
        seqs_not_compat_last_C = ['ATGCATGC', 'CGTACGTC','ATGCATGC', 'CGTACGTC']
        seqs_not_compat_last_G = ['ATGCATGG', 'CGTACGTG','ATGCATGG', 'CGTACGTG']
        seqs_not_compat_last_CG = ['ATGCATGG', 'CGTACGTG','ATGCATGG', 'CGTACGTC']

        self.assertTrue(colour_balanced(seqs_compat_at, channel_2 = False, channel_4 = True))
        self.assertTrue(colour_balanced(seqs_compat_gc, channel_2 = False, channel_4 = True))
        self.assertTrue(colour_balanced(seqs_compat_mix, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_A, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_T, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_AT, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_G, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_C, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_first_GC, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_A, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_T, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_mid_AT, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_A, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_T, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_AT, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_C, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_G, channel_2 = False, channel_4 = True))
        self.assertFalse(colour_balanced(seqs_not_compat_last_CG, channel_2 = False, channel_4 = True))

    def test_colour_balanced_indices(self):
        seqs = enumerate_all_sequences(seqLen = 4, alphabet = ['A', 'T', 'G', 'C'])

        good_4 = generate_colour_balanced_indices(seqs, 4, 10000 , channel_2 = True, channel_4 = True)

        for i in range(0, len(good_4)):
            assert colour_balanced(good_4[i], channel_2 = True, channel_4 = True)
            assert colour_balanced(good_4[i], channel_2 = True, channel_4 = False)
            assert colour_balanced(good_4[i], channel_2 = False, channel_4 = True)


        seqs_bad = enumerate_all_sequences(seqLen = 4, alphabet = ['T', 'G'])
        #Test all combos then fail
        no_good_seqs = generate_colour_balanced_indices(seqs_bad, 4 , 10000 , channel_2 = True, channel_4 = True)
        self.assertEqual(len(no_good_seqs), 0)
        # Test what happens with too few seqs
        too_few_bad_seqs = generate_colour_balanced_indices(seqs_bad, 32 , 10000 , channel_2 = True, channel_4 = True)
        self.assertEqual(len(too_few_bad_seqs), 0)

        too_few_good_seqs = generate_colour_balanced_indices(seqs, 4**4+1 , 10000 , channel_2 = True, channel_4 = True)
        self.assertEqual(len(too_few_good_seqs), 0)

class TestOutput(unittest.TestCase):
    # Grab all indexes from index directories.
    # Ordered by plate number.
    def get_all_seqs_from_dir(self, dirname):
        indexes = []
        for file in sorted(os.listdir(dirname)):
            with open(f"{dirname}{file}", "r") as f:
                lines = f.readlines()
                lines_strip = [line.strip("\n") for line in lines]
                lines_split = [line.split("\t") for line in lines_strip]
                for line in lines_split:
                    for seq in line:
                        indexes.append(seq)
        return indexes

    # Grab final primers from IDT order.
    # Copied these from the final order sheet.
    def get_primers(self, filename):
        with open(filename, "r") as f:
            lines = f.readlines()
            return [line.strip("\n") for line in lines]

    def get_indexes_from_primers(self, primers, regex):
        return [re.match(regex, primer).group(1) for primer in primers]

    # Check that all pairwise distances between sequences in
    # seqs are at least dist.
    def check_distance(self, seqs, dist):
        passed = True
        failed_seqs = []
        for i, seq in enumerate(seqs):
            for comp in seqs[:i] + seqs[i+1:]:
                if levenshtein_dist(seq, comp) < dist:
                    passed = False
                    failed_seqs.append((seq, comp))
        # if not passed:
        #     print(f"Failed! {failed_seqs}")
        return passed

    # Check that each group of 4 index sequences
    # (starting at 0) is colour balanced
    def check_colour_balance(self, seqs):
        balanced = True
        unbalanced = []
        for i in range(0, len(seqs), 4):
            group = seqs[i: i + 4]
            if not colour_balanced(group):
                balanced = False
                unbalanced.append(group)
        # if not balanced:
        #     print("Failed!")
        #     print(unbalanced)
        return balanced

    def lists_equal(self, l1, l2):
        if len(l1) == len(l2) and len(l1) == sum([i == j for i, j in zip(l1, l2)]):
            return True
        else:
            return False

    # tests to make sure that check_distance and check_colour_balance
    # fail when they should fail
    def test_check_distance(self):
        i7s = self.get_all_seqs_from_dir("Final_Indexes/Plates/Indexes/i7/")
        pos = i7s + ["CAGAAGAT"]
        self.assertFalse(self.check_distance(pos, 3))

    def test_check_colour_balance(self):
        i7s = self.get_all_seqs_from_dir("Final_Indexes/Plates/Indexes/i7/")
        i7s[3] = i7s[0]
        i7s[2] = i7s[0]
        self.assertFalse(self.check_colour_balance(i7s))

    def test_generated_indexes(self):
        i7s = self.get_all_seqs_from_dir("Final_Indexes/Plates/Indexes/i7/")
        i5s = self.get_all_seqs_from_dir("Final_Indexes/Plates/Indexes/i5/")
        self.assertTrue(self.check_distance(i7s, 3), msg="All i7 indexes have minimum Lev. distance 3")
        self.assertTrue(self.check_distance(i5s, 3), msg="All i5 indexes have minimum Lev. distance 3")
        self.assertTrue(self.check_colour_balance(i7s), msg="All i7 indexes are colour balanced in groups of 4 across rows")
        self.assertTrue(self.check_colour_balance(i5s), msg="All i7 indexes are colour balanced in groups of 4 across rows")

    def test_primer_sequences(self):
        i7_regex = re.compile("CAAGCAGAAGACGGCATACGAGAT([ATGC]{8})GTGACTGGAGTTCAGACGT\*G")
        i5_regex = re.compile("AATGATACGGCGACCACCGAGATCTACAC([ATGC]{8})ACACTCTTTCCCTACACGAC\*G")
        test_dir = pathlib.Path(__file__).parent
        proj_top_dir = str((test_dir / '..').resolve())
        i7s = self.get_all_seqs_from_dir(f"{proj_top_dir}/Final_Indexes/Plates/Indexes/i7/")
        i5s = self.get_all_seqs_from_dir(f"{proj_top_dir}/Final_Indexes/Plates/Indexes/i5/")
        i7_primers = self.get_primers(f"{proj_top_dir}/test/ordered_primers/i7s.txt")
        i5_primers = self.get_primers(f"{proj_top_dir}/test/ordered_primers/i5s.txt")
        i7_primer_indexes = self.get_indexes_from_primers(i7_primers, i7_regex)
        i5_primer_indexes = self.get_indexes_from_primers(i5_primers, i5_regex)
        self.assertEqual(i7s[:96], i7_primer_indexes, msg="Generated and ordered i7 indexes match")
        self.assertEqual(i5s[:96], i5_primer_indexes, msg="Generated and ordered i5 indexes match")

    # Check that blocklist option works
    # Create a test_list.txt and place at top level directory of project.
    # Run python -m GIL.GIL generate_indexes --blocklist test_list.txt before running the test

    def check_distance_between_sets(self, seqs1, seqs2, dist):
        passed = True
        failed_seqs = []
        for seq1 in seqs1:
            for seq2 in seqs2:
                if levenshtein_dist(seq1, seq2) < dist:
                    passed = False
                    failed_seqs.append((seq1, seq2))
        # if not passed:
        #     print(f"Failed! {failed_seqs}")
        return passed

    # verify that check fails when it should
    def test_check_distance(self):
        self.assertFalse(self.check_distance_between_sets(['ATCGTGAC', 'AAAAAAAA'], ['TCGTGACG', 'AAATAAAA'], 3))
    
    # uncomment this test after running python -m GIL.GIL generate_indexes --blocklist test_list.txt with the blocklist described above
    # def test_blocklist(self):
    #     test_dir = pathlib.Path(__file__).parent
    #     proj_top_dir = str((test_dir / '..').resolve())

    #     with open(f"{proj_top_dir}/test_list.txt", 'r') as f:
    #         blocklist = [index.strip() for index in f.readlines()]

    #     i7s = self.get_all_seqs_from_dir(f"{proj_top_dir}/Output/Plates/Indexes/i7/")
    #     i5s = self.get_all_seqs_from_dir(f"{proj_top_dir}/Output/Plates/Indexes/i5/")

    #     # reverse complement i7s and i5s, as blocklist filters index sequences before they have been reverse complemented
    #     # and placed into primer sequences

    #     i7s = [reverse_complement(seq) for seq in i7s]
    #     i5s = [reverse_complement(seq) for seq in i5s]

    #     self.assertTrue(self.check_distance_between_sets(i7s, blocklist, 3), msg="All i7 indexes have minimum Lev. distance 3 from blocklist")
    #     self.assertTrue(self.check_distance_between_sets(i5s, blocklist, 3), msg="All i5 indexes have minimum Lev. distance 3 from blocklist")

if __name__ == '__main__':
    unittest.main()