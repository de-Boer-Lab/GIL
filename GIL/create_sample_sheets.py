import argparse
import os
import pathlib
import sys
from GIL.tools import *

def make_sample_sheet_UDI(i7s_file, i5s_file, i5_direction, plate_name, dir, template):
  i7s = read_tsv(i7s_file)
  i5s = read_tsv(i5s_file)
  i7s_rc = []
  
  for row in i7s:
    i7s_rc.append([reverse_complement(index) for index in row])
  i5s_rc = []
  for row in i5s:
    i5s_rc.append([reverse_complement(index) for index in row])

  if i5_direction == "forward":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s, plate_name)
  elif i5_direction == "reverse":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s_rc, plate_name)
  else:
    raise ValueError('i5_direction must be "forward" or "reverse"')

  with open(template, 'r') as f:
    sample_template = f.read()
  
  final = sample_template + "\n" + index_lines

  os.makedirs(dir, exist_ok=True)

  with open(f"{dir}/{plate_name}_sample_sheet_{i5_direction}.csv", 'w') as f:
      f.write(final)


def make_sample_sheet_CDI(i7s_file, i5s_file, i7_row, i5_row, i5_direction, plate_name, dir, template):
  
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

  i7s_rc = []
  for row in i7s:
    i7s_rc.append([reverse_complement(index) for index in row])
  i5s_rc = []
  for row in i5s:
    i5s_rc.append([reverse_complement(index) for index in row])

  if i5_direction == "forward":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s, plate_name)
  elif i5_direction == "reverse":
    index_lines = make_sample_sheet_lines(i7s_rc, i5s_rc, plate_name)
  else:
    raise ValueError('i5_direction must be "forward" or "reverse"')

  with open(template, 'r') as f:
    sample_template = f.read()
  
  final = sample_template + "\n" + index_lines

  os.makedirs(dir, exist_ok=True)

  with open(f"{dir}/{plate_name}_sample_sheet_{i5_direction}.csv", 'w') as f:
    f.write(final)


def make_sample_sheet_lines(i7s, i5s, plate_name):
  lines = []
  for i, col in enumerate(zip(zip(*i7s), zip(*i5s))):
    for j, letter in enumerate(["A", "B", "C", "D", "E", "F", "G", "H"]):
      lines.append(f",{plate_name},{letter}{i + 1},{col[0][j]},{col[1][j]}")
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

  GIL_dir = pathlib.Path(__file__).parent
  template_dir = GIL_dir / 'templates' / 'sample_sheet_template.csv'

  if args.unique:
    make_sample_sheet_UDI(args.i7s, args.i5s, i5_direction="forward", dir=args.out_dir, plate_name=args.plate_name, template=template_dir)
    make_sample_sheet_UDI(args.i7s, args.i5s, i5_direction="reverse", dir=args.out_dir, plate_name=args.plate_name, template=template_dir)
  
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
                          i5_direction="forward", dir=args.out_dir, plate_name=args.plate_name, template=template_dir)
    make_sample_sheet_CDI(args.i7s, args.i5s, i7_row=row_num[args.i7_row], i5_row=row_num[args.i5_row],
                          i5_direction="reverse", dir=args.out_dir, plate_name=args.plate_name, template=template_dir)


if __name__ == '__main__':
  main()
