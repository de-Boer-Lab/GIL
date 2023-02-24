import argparse
import os
import pathlib
import sys
from GIL.tools import *


def make_sample_sheet_UDI(i7s_file, i5s_file, plate_name, out_dir):
  i7s = read_tsv(i7s_file)
  i5s = read_tsv(i5s_file)

  info = ('"Forward strand workflow: NovaSeq 6000 with v1.0 reagent kits, '
          'MiSeq, HiSeq 2000/2500, NextSeq 2000 (Sample Sheet v2)"\n'
          '"Reverse strand workflow: iSeq, NovaSeq 6000 with v1.5 reagent kits, '
          'MiniSeq, NextSeq 500/550, HiSeq 3000/4000/X, NextSeq 2000 (Sample Sheet v1)"\n\n')

  header = ("Plate Name,Well,i7 (Read 1) Sequence for Sample Sheet,"
            "i5 (Read 2) Sequence for Sample Sheet (Forward Strand Workflow),"
            "i5 (Read 2) Sequence for Sample Sheet (Reverse Strand Workflow)\n")
  
  index_lines = make_sample_sheet_lines(i7s, i5s, plate_name)

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_sample_sheet.csv", 'w') as f:
    f.write(info + header + index_lines)


def make_sample_sheet_CDI(i7s_file, i5s_file, i7_row, i5_row, plate_name, out_dir):
  
  # Read all indexes in.
  i7s = read_tsv(i7s_file)
  i5s = read_tsv(i5s_file)

  # Keep only the indexes that will
  # be used to make the combinatorial plate.
  i7_subset = i7s[i7_row]
  i5_subset = i5s[i5_row][:8]

  # Copy the i7 subset along cols to
  # make 96 indexes.
  i7s = [i7_subset] * 8

  # Copy the i5 subset along rows to
  # make 96 indexes.
  i5s = [[i] * 12 for i in i5_subset]

  info = ('"Forward strand workflow: NovaSeq 6000 with v1.0 reagent kits, '
          'MiSeq, HiSeq 2000/2500, NextSeq 2000 (Sample Sheet v2)"\n'
          '"Reverse strand workflow: iSeq, NovaSeq 6000 with v1.5 reagent kits, '
          'MiniSeq, NextSeq 500/550, HiSeq 3000/4000/X, NextSeq 2000 (Sample Sheet v1)"\n\n')

  header = ("Plate Name,Well,i7 (Read 1) Sequence for Sample Sheet,"
            "i5 (Read 2) Sequence for Sample Sheet (Forward Strand Workflow),"
            "i5 (Read 2) Sequence for Sample Sheet (Reverse Strand Workflow)\n")
  
  index_lines = make_sample_sheet_lines(i7s, i5s, plate_name)

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_index_sheet.csv", 'w') as f:
    f.write(info + header + index_lines)


def make_sample_sheet_lines(i7s, i5s, plate_name):
  lines = []
  letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
  for letter, i7_row, i5_row in zip(letters, i7s, i5s):
    for i, (i7_sequence, i5_sequence) in enumerate(zip(i7_row, i5_row)):
      lines.append(f"{plate_name},{letter}{i + 1},{reverse_complement(i7_sequence)},{i5_sequence},{reverse_complement(i5_sequence)}")
  return "\n".join(lines)


def main(argv=sys.argv[1:]):
  parser = argparse.ArgumentParser()
  parser.add_argument("--unique", action="store_true", default=False, help="Create sample sheet for unique dual indexes")
  parser.add_argument("--i7s", type=str, help="Path to .tsv containing i7 indexes", required=True)
  parser.add_argument("--i5s", type=str, help="Path to .tsv containing i5 indexes", required=True)
  parser.add_argument("--plate-name", type=str, help="Name of the plate", required=True)
  parser.add_argument("--i7-row", type=str, help="Row that i7 indexes were taken from for creating CDI plate",
                      choices=["A", "B", "C", "D", "E", "F", "G", "H"], required="--unique" not in argv)
  parser.add_argument("--i5-row", type=str, help="Row that i5 indexes were taken from for creating CDI plate",
                      choices=["A", "B", "C", "D", "E", "F", "G", "H"], required="--unique" not in argv)
  parser.add_argument("--out-dir", default="Output/Sample_Sheets/", type=str, help="Path to save sample sheet to")
  args = parser.parse_args(argv)


  if args.unique:
    make_sample_sheet_UDI(args.i7s, args.i5s, dir=args.out_dir, plate_name=args.plate_name)
  
  else:
    row_num = {
      'A': 0,
      'B': 1,
      'C': 2,
      'D': 3,
      'E': 4,
      'F': 5,
      'G': 6,
      'H': 7
    }

    make_sample_sheet_CDI(args.i7s, args.i5s, i7_row=row_num[args.i7_row], i5_row=row_num[args.i5_row],
                          out_dir=args.out_dir, plate_name=args.plate_name)


if __name__ == '__main__':
  main()
