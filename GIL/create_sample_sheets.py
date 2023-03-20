import argparse
import os
import pathlib
import sys
from GIL.tools import *

def get_indexes(index_file):
  try:
    indexes = read_tsv(index_file)
  except ValueError as e:
    print(e)
    print("You likely input an incorrect path to an index TSV file. Check that you provided "
          "correct --i7s and --i5s paths to index TSV files generated by GIL.")
    sys.exit()

  row_lengths = [len(row) for row in indexes]
  if set(row_lengths) != {12} or len(indexes) != 8:
    raise ValueError(f"Error: indexes extracted from {index_file} do not conform to the correct 96 well layout")

  return indexes

def make_index_sheet_lines(i7s, i5s, plate_name):
  lines = []
  letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
  for letter, i7_row, i5_row in zip(letters, i7s, i5s):
    for i, (i7_sequence, i5_sequence) in enumerate(zip(i7_row, i5_row)):
      lines.append(f"{plate_name},{letter}{i + 1},{reverse_complement(i7_sequence)},{i5_sequence},{reverse_complement(i5_sequence)}")
  return "\n".join(lines)


def make_sample_sheet_lines(i7s, i5s, plate_name):
  lines = []
  letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
  for letter, i7_row, i5_row in zip(letters, i7s, i5s):
    for i, (i7_sequence, i5_sequence) in enumerate(zip(i7_row, i5_row)):
      lines.append(f",{plate_name},{letter}{i + 1},{i7_sequence},{i5_sequence}")
  return "\n".join(lines)


def make_index_sheet_UDI(i7s_file, i5s_file, plate_name, out_dir):
  # Read all indexes in.
  try:
    i7s = get_indexes(i7s_file)
    i5s = get_indexes(i5s_file)
  except ValueError as e:
    print(e)
    print("Please ensure that the index TSV files you supplied are formatted in 8 rows by "
          "12 columns of tab-separated indexes.")
    sys.exit()

  info = ('"Forward strand workflow: NextSeq 2000 (Sample Sheet v2), NovaSeq 6000 with v1.0 reagent kits, '
          'MiSeq, HiSeq 2000/2500"\n'
          '"Reverse complement workflow: NextSeq 2000 (Sample Sheet v1), iSeq, NovaSeq 6000 with v1.5 reagent kits, '
          'MiniSeq, NextSeq 500/550, HiSeq 3000/4000/X"\n\n')

  header = ("Plate Name,Well,i7 (Index 1) Sequence for Sample Sheet,"
            "i5 (Index 2) Sequence for Sample Sheet (Forward Strand Workflow),"
            "i5 (Index 2) Sequence for Sample Sheet (Reverse Complement Workflow)\n")
  
  index_lines = make_index_sheet_lines(i7s, i5s, plate_name)

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_index_sheet.csv", 'w') as f:
    f.write(info + header + index_lines)


def make_index_sheet_CDI(i7s_file, i5s_file, i7_row, i5_row, plate_name, out_dir):
  # Read all indexes in.
  try:
    i7s = get_indexes(i7s_file)
    i5s = get_indexes(i5s_file)
  except ValueError as e:
    print(e)
    print("Please ensure that the index TSV files you supplied are formatted in 8 rows by "
          "12 columns of tab-separated indexes.")
    sys.exit()

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

  info = ('"Forward strand workflow: NextSeq 2000 (Sample Sheet v2), NovaSeq 6000 with v1.0 reagent kits, '
          'MiSeq, HiSeq 2000/2500"\n'
          '"Reverse complement workflow: NextSeq 2000 (Sample Sheet v1), iSeq, NovaSeq 6000 with v1.5 reagent kits, '
          'MiniSeq, NextSeq 500/550, HiSeq 3000/4000/X"\n\n')

  header = ("Plate Name,Well,i7 (Index 1) Sequence for Sample Sheet,"
            "i5 (Index 2) Sequence for Sample Sheet (Forward Strand Workflow),"
            "i5 (Index 2) Sequence for Sample Sheet (Reverse Complement Workflow)\n")
  
  index_lines = make_index_sheet_lines(i7s, i5s, plate_name)

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_index_sheet.csv", 'w') as f:
    f.write(info + header + index_lines)


def make_sample_sheet_UDI(i7s_file, i5s_file, i5_direction, plate_name, out_dir, template):
  # Read all indexes in.
  try:
    i7s = get_indexes(i7s_file)
    i5s = get_indexes(i5s_file)
  except ValueError as e:
    print(e)
    print("Please ensure that the index TSV files you supplied are formatted in 8 rows by "
          "12 columns of tab-separated indexes.")
    sys.exit()

  i7s_rc = []
  for row in i7s:
    i7s_rc.append([reverse_complement(index) for index in row])
  i5s_rc = []
  for row in i5s:
    i5s_rc.append([reverse_complement(index) for index in row])

  if i5_direction == "forward_strand":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s, plate_name)
  elif i5_direction == "reverse_complement":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s_rc, plate_name)
  else:
    raise ValueError('i5_direction must be "forward_strand" or "reverse_complement"')

  with open(template, 'r') as f:
    sample_template = f.read()
  
  final = sample_template + "\n" + index_lines

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_sample_sheet_{i5_direction}_workflow.csv", 'w') as f:
      f.write(final)


def make_sample_sheet_CDI(i7s_file, i5s_file, i7_row, i5_row, i5_direction, plate_name,
                          out_dir, template):
  # Read all indexes in.
  try:
    i7s = get_indexes(i7s_file)
    i5s = get_indexes(i5s_file)
  except ValueError as e:
    print(e)
    print("Please ensure that the index TSV files you supplied are formatted in 8 rows by "
          "12 columns of tab-separated indexes.")
    sys.exit()

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

  i7s_rc = []
  for row in i7s:
    i7s_rc.append([reverse_complement(index) for index in row])
  i5s_rc = []
  for row in i5s:
    i5s_rc.append([reverse_complement(index) for index in row])

  if i5_direction == "forward_strand":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s, plate_name)
  elif i5_direction == "reverse_complement":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s_rc, plate_name)
  else:
    raise ValueError('i5_direction must be "forward_strand" or "reverse_complement"')

  with open(template, 'r') as f:
    sample_template = f.read()
  
  final = sample_template + "\n" + index_lines

  os.makedirs(out_dir, exist_ok=True)

  with open(f"{out_dir}/{plate_name}_sample_sheet_{i5_direction}_workflow.csv", 'w') as f:
    f.write(final)


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

  GIL_dir = pathlib.Path(__file__).parent
  template_dir = GIL_dir / 'templates' / 'sample_sheet_template.csv'


  if args.unique:
    make_index_sheet_UDI(args.i7s, args.i5s, out_dir=args.out_dir, plate_name=args.plate_name)
    make_sample_sheet_UDI(args.i7s, args.i5s, i5_direction="forward_strand", out_dir=args.out_dir + "/Forward_Strand_Workflow_Sample_Sheets/", plate_name=args.plate_name, template=template_dir)
    make_sample_sheet_UDI(args.i7s, args.i5s, i5_direction="reverse_complement", out_dir=args.out_dir + "/Reverse_Complement_Workflow_Sample_Sheets/", plate_name=args.plate_name, template=template_dir)
  
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

    make_index_sheet_CDI(args.i7s, args.i5s, i7_row=row_num[args.i7_row], i5_row=row_num[args.i5_row],
                          out_dir=args.out_dir, plate_name=args.plate_name)
    make_sample_sheet_CDI(args.i7s, args.i5s, i7_row=row_num[args.i7_row], i5_row=row_num[args.i5_row], 
                          i5_direction="forward_strand", out_dir=args.out_dir + "/Forward_Strand_Workflow_Sample_Sheets/", plate_name=args.plate_name, template=template_dir)
    make_sample_sheet_CDI(args.i7s, args.i5s, i7_row=row_num[args.i7_row], i5_row=row_num[args.i5_row], 
                          i5_direction="reverse_complement", out_dir=args.out_dir + "/Reverse_Complement_Workflow_Sample_Sheets/", plate_name=args.plate_name, template=template_dir)


if __name__ == '__main__':
  main()
