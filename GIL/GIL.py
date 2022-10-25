'''
GIL.py - Generate indexes for multiplexed sequencing libraries

==============================================================

Choose from two tools:
  - generate_indexes
  - create_sample_sheet

For help with a specific tool, type:
    GIL <tool> --help

To use a specific tool, type:
    GIL <tool> [tool arguments]
'''

import sys
import importlib

def main():
    args = sys.argv

    if len(args) <2 or args[1] == "--help" or args[1] == "-h":
        print("See https://github.com/de-Boer-Lab/IndexingLibrary for full GIL documentation")
        print(globals()["__doc__"])
        return

    command = args[1]

    try:
        module = importlib.import_module("GIL." + command, "GIL")
    except ImportError:
        print(f"{command} is not a GIL command. See GIL -h")
        print("See https://github.com/de-Boer-Lab/IndexingLibrary for full GIL documentation")
        print(globals()["__doc__"])
        return

    module.main(sys.argv[2:])


if __name__ == "__main__":
    sys.exit(main())