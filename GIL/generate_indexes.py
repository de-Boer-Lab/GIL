import argparse
import os
import random
from GIL.create_sample_sheets import *
from GIL.pipeline_components import *
from GIL.tools import *

try:
 import Levenshtein
 levenshtein_dist = Levenshtein.distance
except ModuleNotFoundError:
    print("Levenshtein module (https://github.com/maxbachmann/Levenshtein) not found. "
            "A much slower implementation of Levenshtein distance will be used.")
    levenshtein_dist = levenshtein

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", default=8, type=int, help="length of indexes")
    parser.add_argument("--dist", default=3, type=int, help="minimum Levenshtein distance between indexes")
    parser.add_argument("--min-GC", default=25, type=int, help="GC content of indexes will be strictly higher than this value")
    parser.add_argument("--max-GC", default=75, type=int, help="GC content of indexes will be strictly lower than this value")
    parser.add_argument("--max-homopolymer", default=2, type=int, help="maximum length of homopolymer repeats in indexes")
    parser.add_argument("--max-dinu", default=2, type=int, help="maximum number of dinucleotide repeats in indexes")
    parser.add_argument("--allow-start-G", action="store_true", help="allow indexes to start with G")
    parser.add_argument("--no-filter-self-priming", action="store_true", default=False, help="don't filter indexes for self-priming potential")
    parser.add_argument("--sample-n", type=int, help="number of random indexes to generate at the start")
    parser.add_argument("--library-type", default="TruSeq", choices=["TruSeq", "Nextera", "custom"], type=str, help="library type to use for primer sequences")
    parser.add_argument("--primer-sequences", type=str, help="sequences of custom primers. Usage: --primer_sequences \"<forward 5'> <forward 3'> <reverse 5'> <reverse 3'>\"")
    parser.add_argument("--library-name", type=str, help="library type name that will be included in file names; will replace \"custom\" for custom primer sequences")
    parser.add_argument("--blocklist", type=str, help="path to file containing blocklist indexes")
    parser.add_argument("--no-mod", action="store_true", default=False, help="don't add a '*' (phosphorothioate modification) to primers in order sheet")
    parser.add_argument("--out-dir", default="Output", type=str, help="path to location of program outputs")
    parser.add_argument("--seed", type=int, help="seed for reproducible index generation")
    parser.add_argument("--company", default="IDT", type=str, choices=["IDT", "Thermo", "Sigma"], help="format to use for oligo order sheet")
    
    args = parser.parse_args(argv)

    # Check for correct argument values

    if args.length < 4:
        sys.exit("Error: --length must be at least 4 to generate enough indexes.")
    
    if args.length > 45:
        print(f"Warning: are you sure you need indexes of length {args.length}? "
               "Colour balancing of long indexes (>45) is difficult, resulting in few or no "
               "generated indexes. You can try relaxing filtering parameters or increasing "
               "--sample-n if you really need indexes this long.")
        
    if args.dist < 1:
        sys.exit("Error: --dist must be > 0 to ensure no duplicate indexes.")
    
    if args.dist > args.length:
        sys.exit(f"Error: --dist cannot be greater than --length. --dist is set to {args.dist}, "
                         f"but --length is set to {args.length}.")
        
    if args.dist >= args.length / 2:
        print(f"Warning: are you sure you need a distance of {args.dist} for indexes of length {args.length}? "
               "It is difficult to generate sufficient indexes when distance >= length/2. "
               "You can try relaxing filtering parameters, or increasing "
               "--sample-n (if your index length is greater than 9) if "
               "you really need a distance this large.")
        
    if args.min_GC < 0 or args.min_GC > 100 or args.max_GC < 0 or args.max_GC > 100:
        sys.exit("Error: --min-GC and max-GC must be between 0 and 100.")
        
    if args.min_GC > args.max_GC:
        sys.exit("Error: you have requested a minimum GC content which is higher "
                 "than the requested maximum GC content. Please ensure --min-GC < --max-GC.")
        
    if args.max_homopolymer < 1 or args.max_homopolymer > args.length // 2:
        sys.exit("Error: --max-homopolymer must be > 1 and < --length // 2.")

    if args.max_dinu < 1 or args.max_dinu > args.length // 2:
        sys.exit("Error: --max-dinu must be > 1 and < --length // 2.")

    if args.seed:
        random.seed(args.seed)

    # Set up primer sequences
    if args.library_type == "TruSeq":
        primer_i5_start = "AATGATACGGCGACCACCGAGATCTACAC"
        primer_i5_end = "ACACTCTTTCCCTACACGACG"
        primer_i7_start = "CAAGCAGAAGACGGCATACGAGAT"
        primer_i7_end = "GTGACTGGAGTTCAGACGTG"
    elif args.library_type == "Nextera":
        primer_i5_start = "AATGATACGGCGACCACCGAGATCTACAC"
        primer_i5_end = "TCGTCGGCAGCGTC"
        primer_i7_start = "CAAGCAGAAGACGGCATACGAGAT"
        primer_i7_end = "GTCTCGTGGGCTCGG"
    else:
        if args.primer_sequences is None:
            sys.exit("Error: primer sequences must be supplied with the --primer_sequences argument "
                            "if custom library type is chosen")
        primer_seqs = args.primer_sequences.split()
        primer_i5_start = primer_seqs[0]
        primer_i5_end = primer_seqs[1]
        primer_i7_start = primer_seqs[2]
        primer_i7_end = primer_seqs[3]

    # Set up library name
    if args.library_name is None:
        library_name = args.library_type
    else:
        library_name = args.library_name

    # Set up blocklist
    if args.blocklist is not None:
        with open(args.blocklist, 'r') as f:
            blocklist = [index.strip() for index in f.readlines()]

    # Generate all possible indexes. For index lengths > 9, generating all possible sequences and filtering them
    # all takes too long. Instead, generate a sample of 5k sequences.
    if args.length < 10 and args.sample_n is None:
        seqs = enumerate_all_sequences(["A", "C", "G", "T"], args.length)
    else:
        if args.sample_n is None:
            sample_n = 5000
        else:
            sample_n = args.sample_n
        GC_weight = (args.max_GC + args.min_GC) / 200
        AT_weight = 1 - GC_weight
        seqs = ["".join(random.choices(["A", "C", "G", "T"], k=args.length, weights=[AT_weight, GC_weight, GC_weight, AT_weight])) for i in range(sample_n)]

    # Filter generated sequences
    if args.allow_start_G:
        seqs = filter_startingG(seqs, 2)
    else:
        seqs = filter_startingG(seqs, 1)
    seqs = filter_GC(seqs, minGC=args.min_GC, maxGC=args.max_GC)
    seqs = filter_homopolymer(seqs, max_repeat=args.max_homopolymer)
    seqs = filter_dinucleotide_repeats(seqs, max_dinu=args.max_dinu)
    if args.blocklist is not None:
        seqs = filter_blocklist(seqs, blocklist, dist=args.dist, dist_fun=levenshtein_dist)
    if not args.no_filter_self_priming:
        seqs_i7 = filter_self_priming(seqs, min_dist=3, primer_5prime=primer_i7_start, primer_3prime=primer_i7_end)
        seqs_i5 = filter_self_priming(seqs, min_dist=3, primer_5prime=primer_i5_start, primer_3prime=primer_i5_end)
    else:
        seqs_i7 = seqs
        seqs_i5 = seqs
    if args.allow_start_G:
        seqs_i5 = filter_endingC(seqs_i5, 2)
    else:
        seqs_i5 = filter_endingC(seqs_i5, 1)

    # Select dissimilar subset and group for colour balance
    seqs_i7 = select_dissimilar_sequences(seqs_i7, dist=args.dist, dist_fun=levenshtein_dist)
    seqs_i5 = select_dissimilar_sequences(seqs_i5, dist=args.dist, dist_fun=levenshtein_dist)
    final_i7 = generate_colour_balanced_indices(seqs_i7, 4, 1_000_000)
    final_i5 = generate_colour_balanced_indices(seqs_i5, 4, 1_000_000)

    # Create index and primer files and order and sample sheets
    try:
        create_primer_plates(final_i7, final_i5, i5_start=primer_i5_start, i5_end=primer_i5_end,
                             i7_start=primer_i7_start, i7_end=primer_i7_end, dir=args.out_dir, plate_name=library_name)
    except ValueError:
        sys.exit("Not enough indexes to create a 96 well plate. You can try again with relaxed parameters, or try "
                 "increasing --sample-n if the length of your indexes is greater than 9.")

    primer_dir = f"{args.out_dir}/Plates/Primers/"
    for plate_file in os.listdir(primer_dir):
        plate_name = plate_file.split(".")[0]
        plate_name_split = plate_name.split("_")
        primer_name = "_".join([plate_name_split[i] for i in [3, 4, 0, 1]])
        make_order_sheet(f"{primer_dir}{plate_file}", plate_name, primer_name, company=args.company, dir=args.out_dir, mod=not args.no_mod)

    GIL_dir = pathlib.Path(__file__).parent
    template_dir = GIL_dir / 'templates' / 'sample_sheet_template.csv'

    index_dir = f"{args.out_dir}/Plates/Indexes/"
    i7_files = sorted(os.listdir(f"{index_dir}i7"))
    i5_files = sorted(os.listdir(f"{index_dir}i5"))
    for i7_file, i5_file in zip(i7_files, i5_files):
        plate_name = i7_file.split(".")[0]
        split_name = plate_name.split("_")
        plate_name = "_".join([split_name[i] for i in [0, 3, 4]])
        make_index_sheet_UDI(f"{index_dir}i7/{i7_file}", f"{index_dir}i5/{i5_file}",
                            plate_name=plate_name, out_dir=f"{args.out_dir}/Sample_Sheets/")
        make_sample_sheet_UDI(f"{index_dir}i7/{i7_file}", f"{index_dir}i5/{i5_file}", plate_name=plate_name, i5_direction="forward_strand",
                              out_dir=f"{args.out_dir}/Sample_Sheets/Forward_Strand_Workflow_Sample_Sheets/", template=template_dir)
        make_sample_sheet_UDI(f"{index_dir}i7/{i7_file}", f"{index_dir}i5/{i5_file}", plate_name=plate_name, i5_direction="reverse_complement",
                              out_dir=f"{args.out_dir}/Sample_Sheets/Reverse_Complement_Workflow_Sample_Sheets/", template=template_dir)

if __name__ == '__main__':
    main()