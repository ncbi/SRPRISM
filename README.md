# SRPRISM - Single Read Paired Read Indel Substitution Minimizer
Version 3.3.3-pre

For questions regarding SRPRISM, please contact
    Aleksandr Morgulis (morgulis@ncbi.nlm.nih.gov)
    or
    Richa Agarwala (agarwala@ncbi.nlm.nih.gov)

## Linux Compilation

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

## Windows Compilation

    Windows 64-bit version of SRPRISM was tested with MS Visual Studio 2022 and
    MS Visual Studio 2022 Community Edition. In order to compile the following
    is required:

        - MS Visual Studio 2022 (or Community Edition of it);
        - if cloning from github: git with git bash (available here: https://git-scm.com/download/win)

    To compile:

        - select a work directory (in the following $SRP_HOME is used as a
            full pathname of that directory);
        - to clone the latest version: start git bash and issue the following commands:
            cd $SRP_HOME
            git clone https://github.com/ncbi/SRPRISM
        - to use the latest stable release: download the latest .zip archive from
          github.com/ncbi/SRPRISM/releases and unpack it in $SRP_HOME
        - start MS Visual Studio 2022 and open solution: $SRP_HOME\SRPRISM\windows\SRPRISM\SRPRISM.sln
        - select Release/x64 configuration
        - in the "Solution Explorer" pane, right click on the solution and select "build solution" in
          the context menu
        - the final executable SRPRISM.exe will be in $SRP_HOME\SRPRISM\windows\SRPRISM\x64\Release\

## Usage

    Please see srprism/README file for usage infromation and examples.

