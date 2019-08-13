# SRPRISM - Single Read Paired Read Indel Substitution Minimizer
Version 3.0.1

For questions regarding SRPRISM, please contact
    Aleksandr Morgulis (morgulis@ncbi.nlm.nih.gov)
    or
    Richa Agarwala (agarwala@ncbi.nlm.nih.gov)

## Compilation

    Download current source code for SRPRISM
    $ git clone https://github.com/ncbi/SRPRISM

    Do following:
    $ cd SRPRISM/srprism

    To enable support for NGS library, needed to access NCBI SRA archive directly,
    please uncomment the following line in the beginning of the top level Makefile
    # export USE_SRA = 1
    This will trigger automatic download and compilation of the NGS library.
    In order to use pre-installed NGS and VDB libraries, uncomment the following
    lines in the top level Makefile and change them to point to the relevant paths
    # export NGS_PATH =
    # export VDB_PATH =

    To build SRPRISM application, do
    $ make

    After successful build, srprism executable can be found in app/ subdirectory.

## Usage

    Please see srprism/README file for usage infromation and examples.

