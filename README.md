# GIL: Generate Indexes for Libraries

GIL is a Python package that designs indexing primers for multiplexed sequencing.

See the [wiki](https://github.com/de-Boer-Lab/GIL/wiki) for explanations of what GIL does and
examples of how to use it. Can't find an answer to your question in the wiki? Start a discussion on GIL's [discussion page](https://github.com/de-Boer-Lab/GIL/discussions).

A web app version of GIL is available here: https://dbl-gil.streamlitapp.com/. If you encounter a message saying  the app is down/booting, please be patient and try again in 30 minutes. This can happen due to inactivitiy. 

![](images/novaseq_gil@0.25x.png)
# Installation

GIL requires Python 3.6 or higher and is not compatible with Python 2. The optional
[Levenshtein](https://pypi.org/project/Levenshtein/) module (`pip install Levenshtein`) makes index generation much faster but is not required.

## Installing from GitHub

To install GIL from GitHub, use

`pip install git+https://github.com/de-Boer-Lab/GIL`

GIL can then be run as follows:

`GIL [COMMAND] [ARGUMENTS]`

For help:

`GIL -h`

or  

`GIL [COMMAND] -h`

See [usage examples](#usage-examples) below for more details, and the [command line arguments](#command-line-arguments) section for explanations of all arguments.

## Running from a cloned repo

If you'd prefer not to install GIL, you can instead clone the repo:

`git clone https://github.com/de-Boer-Lab/GIL.git`

Commands can only be run from the top-level GIL directory as follows:

`python -m GIL.[COMMAND] [ARGUMENTS]`

For example, to generate indexes of length 10, run

`python -m GIL.generate_indexes --length 10`

#