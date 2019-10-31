/*  $Id: srprism.cpp 590234 2019-07-25 16:28:25Z morgulis $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Aleksandr Morgulis, Denis Vakatov, Vladimir Ivanov
 *
 * File Description: srpsism - short read aligner
 *
 */

#include <common/def.h>

#include <iostream>
#include <algorithm>
#include <memory>

#include <common/exception.hpp>
#include <options_parser.hpp>
#include <common/trace.hpp>
#include <common/util.hpp>

#include <srprism/search.hpp>
#include <srprism/mkidx.hpp>
#include <srprism/out_tabular.hpp>
#include <srprism/out_sam.hpp>

USE_NS( STD_SCOPES::common )
USE_NS( STD_SCOPES::srprism )

//------------------------------------------------------------------------------
// program info
//
const char * PROGRAM_DESCRIPTION = 
"\
Fast Short Read Aligner \
";

//------------------------------------------------------------------------------
// positional parameters info
//
static const COptionsParser::TParamPos CMD_POS = 0;
static const COptionsParser::TParamPos HLP_SEC_POS = 1;
static const std::string CMD_LABEL = "cmd";
static const std::string CMD_DESCR = R"(
    Action to perform. Possible values are:
          help [search|mkindex] - get usage help;
                                  general help if no option is given;
                                  otherwise help on specified command.
          search                - search for occurrences of the
                                  queries in the database;
          mkindex               - create index from a source database.
    Type 'srprism help search' for more help on search command.
    Type 'srprism help mkindex' for more help on mkindex comamnd.
)";

/*
static const std::string CMD_DESCR = "\
\tAction to perform. Possible values are:\n\
\t\t* help [search|mkindex] - \tget usage help;\n\
\t\t                          \t\tgeneral help if no option is given;\
\t\t                          \t\totherwise help on specified command.\n\
\t\t* search                - \tsearch for occurrences of the \
\t\t                          \tqueries in the database;\n\
\t\t* mkindex               - \tcreate index from a source database.\n\
\tType 'srprism help search' for more help on search command.\n\
\tType 'srprism help mkindex' for more help on mkindex comamnd.\n\
";
*/

//------------------------------------------------------------------------------
// common options
//
static const std::string TRACE_LEVEL_KEY     = "trace-level";
static const std::string TRACE_LEVEL_SKEY    = "";
static const std::string TRACE_LEVEL_LABEL   = "min_level";
static const std::string TRACE_LEVEL_DEFAULT = "warning";
static const std::string TRACE_LEVEL_DESCR   = "\
\tMinimum message level to report to the log stream. Possible values are \
\"debug\", \"info\", \"warning\", \"error\", \"quiet\".\n\
";

static const std::string LOG_FNAME_KEY   = "log-file";
static const std::string LOG_FNAME_SKEY  = "";
static const std::string LOG_FNAME_LABEL = "file-name";
static const std::string LOG_FNAME_DESCR = "\
\tFile for storing diagnostic messages. Default is standard error.\n\
";

//------------------------------------------------------------------------------
// search options
//
static const std::string SEARCH_INDEX_KEY   = "index";
static const std::string SEARCH_INDEX_SKEY  = "I";
static const std::string SEARCH_INDEX_LABEL = "basename";
static const std::string SEARCH_INDEX_DESCR = "\
\tBase name for database index files.\n\
";

static const std::string SEARCH_INPUT_KEY   = "input";
static const std::string SEARCH_INPUT_SKEY  = "i";
static const std::string SEARCH_INPUT_LABEL = "input-file(s)";
static const std::string SEARCH_INPUT_DESCR = "\
\tSpecifies the source of the queries. The exact format depends on the \
value of \"--input-format\" option.\n\
";

static const std::string SEARCH_OUTPUT_KEY   = "output";
static const std::string SEARCH_OUTPUT_SKEY  = "o";
static const std::string SEARCH_OUTPUT_LABEL = "output-stream";
static const std::string SEARCH_OUTPUT_DESCR = "\
\tFile name for program output. By default the standard output stream will be \
used.\n\
";

static const std::string SEARCH_MODE_KEY     = "mode";
static const std::string SEARCH_MODE_SKEY    = "m";
static const std::string SEARCH_MODE_LABEL   = "search-mode";
static const std::string SEARCH_MODE_DEFAULT = "min-err";
static const std::string SEARCH_MODE_DESCR   = "\
\tSearch mode; possible values are \"min-err\", \"bound-err\", \
\"sum-err\", \"partial\".\n\
";

static const std::string SEARCH_NERR_KEY     = "errors";
static const std::string SEARCH_NERR_SKEY    = "n";
static const std::string SEARCH_NERR_LABEL   = "max-errors";
static const std::string SEARCH_NERR_DEFAULT = "5";
static const std::string SEARCH_NERR_DESCR   = "\
\tSearch for alignments with at most this many errors (maximum is 15).\n\
";

static const std::string SEARCH_NRES_KEY     = "results";
static const std::string SEARCH_NRES_SKEY    = "r";
static const std::string SEARCH_NRES_LABEL   = "max-results";
static const std::string SEARCH_NRES_DEFAULT = "10";
static const std::string SEARCH_NRES_DESCR   = "\
\tMaximum number of results to report per query (maximum is 255).\n\
";

static const std::string SEARCH_DIST_KEY     = "pair-distance";
static const std::string SEARCH_DIST_SKEY    = "s";
static const std::string SEARCH_DIST_LABEL   = "target-pair-distance";
static const std::string SEARCH_DIST_DEFAULT = "500";
static const std::string SEARCH_DIST_DESCR   = "\
\tFor paired search, the target distance between pair alignments in bases.\n\
";

static const std::string SEARCH_FUZZ_KEY     = "pair-distance-fuzz";
static const std::string SEARCH_FUZZ_SKEY    = "f";
static const std::string SEARCH_FUZZ_LABEL   = "fuzz-value";
static const std::string SEARCH_FUZZ_DEFAULT = "490";
static const std::string SEARCH_FUZZ_DESCR   = "\
\tMaximum acceptable radius (in bases) around the target distance between \
paired alignments. If d is the value of \"--pair-distance\" option and f is \
the value of the \"--pair-distance-fuzz\" option than the distance between \
the alignments in a pair should be in the interval [d-f,d+f] for the \
pair to be reported.\n\
";

static const std::string SEARCH_INFMT_KEY     = "input-format";
static const std::string SEARCH_INFMT_SKEY    = "F";
static const std::string SEARCH_INFMT_LABEL   = "format-name";
static const std::string SEARCH_INFMT_DEFAULT = "fasta";
static const std::string SEARCH_INFMT_DESCR   = "\
\tThe input format name. The possible values are \
\"fasta\", \"fastq\", \"cfasta\", \"cfastq\", \"sam\", \"sra\". \
See the software documentation for the details of different supported \
input formats.\n\
\tIn the case of paired queries some formats use two files that are \
read in parallel by srprizm. In this case reading from standard input \
is impossible. Also, in this case, the value of \"--input <name>\" \
option is interpreted as follows:\n\
\t\t* if name is of the form \"<name1>,<name2>\" then name1 and name2 \
are assumed to be the files containing respectivel first and second \
members of the pairs.\n\
\tBoth files have to contain the same number of sequences and the ids of \
of the corresponding sequences should be identical.\n\
\tThe formats that have the above restriction are \"fasta\", \"fastq\", \"cfasta\", \"cfastq\".\n\
\tFor formats \"fasta\", \"fastq\", \"cfasta\", \"cfastq\", if \
\ta paired search is requested but only a single input file is specified, \
\tthen that file has to contain even number of entries, with each consequtive \
\tpair considered one paired-end query.\n\
";

static const std::string SEARCH_ICOMPR_KEY  = "input-compression";
static const std::string SEARCH_ICOMPR_SKEY = "";
static const std::string SEARCH_ICOMPR_LABEL = "compression-type";
static const std::string SEARCH_ICOMPR_DEFAULT = "auto";
static const std::string SEARCH_ICOMPR_DESCR   = "\
\tCompression type used for input. The possible values are \
\"auto\" (default), \"none\", \"gzip\", and \"bzip2\". If the value \
given is \"auto\" then the type of compression is guessed from the file \
extension.\n\
";

/*
static const std::string SEARCH_OUTFMT_KEY     = "output-format";
static const std::string SEARCH_OUTFMT_SKEY    = "O";
static const std::string SEARCH_OUTFMT_LABEL   = "format-name";
static const std::string SEARCH_OUTFMT_DEFAULT = "tabular";
static const std::string SEARCH_OUTFMT_DESCR   = "\
\tThe output format name. The possible values are \"tabular\", \"sam\". \
See the software documentation for the details of different supported \
output formats.\n\
";
*/

static const std::string SEARCH_MEM_KEY     = "memory";
static const std::string SEARCH_MEM_SKEY    = "M";
static const std::string SEARCH_MEM_LABEL   = "megabytes";
static const std::string SEARCH_MEM_DEFAULT = "4096";
static const std::string SEARCH_MEM_DESCR   = "\
\tDo not use more than this many megabytes of memory for internal \
dynamic data structures. This number does not include the footprint \
of the executable code, static data, or stack.\n\
";

static const std::string SEARCH_BATCH_KEY     = "batch";
static const std::string SEARCH_BATCH_SKEY    = "b";
static const std::string SEARCH_BATCH_LABEL   = "batch_size";
static const std::string SEARCH_BATCH_DESCR   = "\
\tProcess input in batches of at most this many queries.\n\
";

static const std::string SEARCH_SBATCH_KEY     = "batch-start";
static const std::string SEARCH_SBATCH_SKEY    = "";
static const std::string SEARCH_SBATCH_LABEL   = "first batch";
static const std::string SEARCH_SBATCH_DESCR   = "\
\tFirst batch to process.\n\
";

static const std::string SEARCH_EBATCH_KEY     = "batch-end";
static const std::string SEARCH_EBATCH_SKEY    = "";
static const std::string SEARCH_EBATCH_LABEL   = "last batch";
static const std::string SEARCH_EBATCH_DESCR   = "\
\tLast batch to process.\n\
";

static const std::string SEARCH_QID_KEY   = "no-qids";
static const std::string SEARCH_QID_SKEY  = "D";
static const std::string SEARCH_QID_DESCR = "\
\tDo not report ids for queries. Use their ordinal number instead.\n\
";

static const std::string SEARCH_SID_KEY   = "no-sids";
static const std::string SEARCH_SID_SKEY  = "";
static const std::string SEARCH_SID_DESCR = "\
\tDo not report ids for database sequences. Use their ordinal number instead.\n\
";

static const std::string SEARCH_TMPDIR_KEY     = "tmpdir";
static const std::string SEARCH_TMPDIR_SKEY    = "T";
static const std::string SEARCH_TMPDIR_LABEL   = "dir-name";
static const std::string SEARCH_TMPDIR_DEFAULT = ".";
static const std::string SEARCH_TMPDIR_DESCR   = "\
\tDirectory to store temporary files.\n\
";

static const std::string SEARCH_PAIRED_LOG_KEY  = "plog";
static const std::string SEARCH_PAIRED_LOG_SKEY = "";
static const std::string SEARCH_PAIRED_LOG_LABEL = "file-name";
static const std::string SEARCH_PAIRED_LOG_DESCR = "\
\tFile name where the ordinal ids of queries participating in \
paired search passes are written.\n\
";

static const std::string SEARCH_PAIRED_KEY   = "paired";
static const std::string SEARCH_PAIRED_SKEY  = "p";
static const std::string SEARCH_PAIRED_LABEL = "search-type";
static const std::string SEARCH_PAIRED_DESCR = "\
\tIf \"true\", force paired search; if \"false\", force unpaired search.\n\
";

static const std::string SEARCH_SKIP_UNMAPPED_KEY = "skip-unmapped";
static const std::string SEARCH_SKIP_UNMAPPED_SKEY = "S";
static const std::string SEARCH_SKIP_UNMAPPED_LABEL = "true|false";
static const std::string SEARCH_SKIP_UNMAPPED_DEFAULT = "true";
static const std::string SEARCH_SKIP_UNMAPPED_DESCR = "\
\tIf \"true\", do not generate records for unmapped queries in SAM output.\n\
";

static const std::string SEARCH_REPEAT_KEY     = "repeat-threshold";
static const std::string SEARCH_REPEAT_SKEY    = "R";
static const std::string SEARCH_REPEAT_LABEL   = "intval";
static const std::string SEARCH_REPEAT_DEFAULT = "4096";
static const std::string SEARCH_REPEAT_DESCR   = "\
\tProcess queries that start with 16-mers that appear at most this many \
times in the database.\n\
";

static const std::string SEARCH_RESCONF_KEY     = "result-conf";
static const std::string SEARCH_RESCONF_SKEY    = "c";
static const std::string SEARCH_RESCONF_LABEL   = "result_configuration_spec";
static const std::string SEARCH_RESCONF_DEFAULT = "0100";
static const std::string SEARCH_RESCONF_DESCR   = "\
\tSelect which paired result configurations should be reported.\n\
\tSelection is specified as a sequence of 4 digits from {0,1}, where \
i-th digit is 1 if the i-th configuration is selected. Configurations \
are encoded by the second mate direction and position in the result \
relative to the first mate, assuming that the first mate is matched \
in forward direction.\n\
\tConfigurations are encoded via the following table:\n\
\t\t0 - the second mate is forward aligned to the right of the first mate;\n\
\t\t1 - the second mate is reverse aligned to the right of the first mate;\n\
\t\t2 - the second mate is forward aligned to the left of the first mate;\n\
\t\t3 - the second mate is reverse aligned to the left of the first mate.\n\
\tThe following aliases are supported:\n\
\t\tillumina = 0100\n\
\t\t454      = 0100\n\
\t\tsolid    = 0010\n\
";

static const std::string SEARCH_SA_START_KEY     = "sa-start";
static const std::string SEARCH_SA_START_SKEY    = "";
static const std::string SEARCH_SA_START_LABEL   = "seeding_area_start_pos";
static const std::string SEARCH_SA_START_DEFAULT = "1";
static const std::string SEARCH_SA_START_DESCR   = "\
\tMake sure that the seeding area is selected so that it starts at or after \
this query position.\n\
";

static const std::string SEARCH_SA_END_KEY     = "sa-end";
static const std::string SEARCH_SA_END_SKEY    = "";
static const std::string SEARCH_SA_END_LABEL   = "seeding_area_end_pos";
static const std::string SEARCH_SA_END_DEFAULT = "8096";
static const std::string SEARCH_SA_END_DESCR   = "\
\tMake sure that the seeding area is selected so that it ends at or before \
this query position.\n\
";

static const std::string SEARCH_EXTRA_TAGS_KEY     = "extra-tags";
static const std::string SEARCH_EXTRA_TAGS_SKEY    = "";
static const std::string SEARCH_EXTRA_TAGS_LABEL   = "string";
static const std::string SEARCH_EXTRA_TAGS_DEFAULT = "";
static const std::string SEARCH_EXTRA_TAGS_DESCR   = "\
\tString (normally fixed extra tags) to add to each SAM output record.\n\
";

static const std::string SEARCH_SD_KEY   = "discover-insert";
static const std::string SEARCH_SD_SKEY  = "";
static const std::string SEARCH_SD_DESCR = "\
\tTry automatically discover the separation between mates in paired-end runs.\n\
";

static const std::string SEARCH_SD_STOP_KEY = "discover-insert-and-stop";
static const std::string SEARCH_SD_STOP_SKEY = "";
static const std::string SEARCH_SD_STOP_DESCR = "\
\tUnconditionally stop after insert size analysis, dump the \
histogram and generate single results.\n\
";

static const std::string SEARCH_SD_HNAME_KEY     = "discover-insert-hist-name";
static const std::string SEARCH_SD_HNAME_SKEY    = "";
static const std::string SEARCH_SD_HNAME_LABEL   = "file-name";
static const std::string SEARCH_SD_HNAME_DEFAULT = "hist.out";
static const std::string SEARCH_SD_HNAME_DESCR = "\
\tFile name to dump the histogram to after the insert size analysis.\n\
";

static const std::string SEARCH_RANDOMIZE_KEY     = "randomize";
static const std::string SEARCH_RANDOMIZE_SKEY    = "";
static const std::string SEARCH_RANDOMIZE_LABEL   = "true|false";
static const std::string SEARCH_RANDOMIZE_DEFAULT = "false";
static const std::string SEARCH_RANDOMIZE_DESCR   = "\
\tTry to randomize the results of each query on subject coordinate.\n\
";

static const std::string SEARCH_RANDOM_SEED_KEY     = "random-seed";
static const std::string SEARCH_RANDOM_SEED_SKEY    = "";
static const std::string SEARCH_RANDOM_SEED_LABEL   = "true|false";
static const std::string SEARCH_RANDOM_SEED_DEFAULT = "false";
static const std::string SEARCH_RANDOM_SEED_DESCR   = "\
\tUse random seed for results randomization (warning: results \
will not be reproducible).\n\
";

static const std::string SEARCH_SAM_HEADER_KEY      = "sam-header";
static const std::string SEARCH_SAM_HEADER_SKEY     = "";
static const std::string SEARCH_SAM_HEADER_LABEL    = "true|false";
static const std::string SEARCH_SAM_HEADER_DEFAULT  = "false";
static const std::string SEARCH_SAM_HEADER_DESCR    = "\
\tPrint the standard SAM header.\n\
";

#ifndef NDEBUG

static const std::string SEARCH_FIX_HC_KEY = "hc";
static const std::string SEARCH_FIX_HC_SKEY = "";
static const std::string SEARCH_FIX_HC_LABEL = "hc_value";
static const std::string SEARCH_FIX_HC_DESCR = "\
\tUse fixed hash configuration for all queries \
(this option if for debugging/testing purposes only: use with caution).\n\
";

#endif

//------------------------------------------------------------------------------
// mkindex options
//
static const std::string MKINDEX_INPUT_KEY   = "input";
static const std::string MKINDEX_INPUT_SKEY  = "i";
static const std::string MKINDEX_INPUT_LABEL = "input-file(s)";
static const std::string MKINDEX_INPUT_DESCR = "\
\tSource database (may be a comma-separated list of file names). \
This options takes precedence over --input-list.\n\
";

static const std::string MKINDEX_INPUT_LIST_KEY = "input-list";
static const std::string MKINDEX_INPUT_LIST_SKEY = "l";
static const std::string MKINDEX_INPUT_LIST_LABEL = "file-name";
static const std::string MKINDEX_INPUT_LIST_DESCR = "\
\tName of the file containing a list of input file names, one \
name per line.\n\
";

static const std::string MKINDEX_ALSPEC_KEY = "alt-loc";
static const std::string MKINDEX_ALSPEC_SKEY = "a";
static const std::string MKINDEX_ALSPEC_LABEL = "file-name";
static const std::string MKINDEX_ALSPEC_DESCR = "\
\tName of the alternative loci specification file.\n\
";

static const std::string MKINDEX_OUTPUT_KEY   = "output";
static const std::string MKINDEX_OUTPUT_SKEY  = "o";
static const std::string MKINDEX_OUTPUT_LABEL = "output-stream";
static const std::string MKINDEX_OUTPUT_DESCR = "\
\tBase name for generated database index files.\n\
";

static const std::string MKINDEX_INFMT_KEY     = "input-format";
static const std::string MKINDEX_INFMT_SKEY    = "F";
static const std::string MKINDEX_INFMT_LABEL   = "format-name";
static const std::string MKINDEX_INFMT_DEFAULT = "fasta";
static const std::string MKINDEX_INFMT_DESCR   = "\
\tThe input format name. The possible values are \"fasta\", \"fastq\", \"cfasta\", \"cfastq\".\n\
";

static const std::string MKINDEX_ICOMPR_KEY  = "input-compression";
static const std::string MKINDEX_ICOMPR_SKEY = "";
static const std::string MKINDEX_ICOMPR_LABEL = "compression-type";
static const std::string MKINDEX_ICOMPR_DEFAULT = "auto";
static const std::string MKINDEX_ICOMPR_DESCR   = "\
\tCompression type used for input. The possible values are \
\"auto\" (default), \"none\", \"gzip\", and \"bzip2\". If the value \
given is \"auto\" then the type of compression is guessed from the file \
extension.\n\
";

static const std::string MKINDEX_OUTFMT_KEY     = "output-format";
static const std::string MKINDEX_OUTFMT_SKEY    = "O";
static const std::string MKINDEX_OUTFMT_LABEL   = "format-name";
static const std::string MKINDEX_OUTFMT_DEFAULT = "standard";
static const std::string MKINDEX_OUTFMT_DESCR   = "\
\tThe output format name. The possible values are \"standard\".\n\
";

static const std::string MKINDEX_MEM_KEY     = "memory";
static const std::string MKINDEX_MEM_SKEY    = "M";
static const std::string MKINDEX_MEM_LABEL   = "megabytes";
static const std::string MKINDEX_MEM_DEFAULT = "4096";
static const std::string MKINDEX_MEM_DESCR   = "\
\tDo not use more than this many megabytes of memory for internal \
dynamic data structures. This number does not include the footprint \
of the executable code, static data, or stack.\n\
";

static const std::string MKINDEX_SEGLEN_KEY     = "seg-letters";
static const std::string MKINDEX_SEGLEN_SKEY    = "";
static const std::string MKINDEX_SEGLEN_LABEL   = "bases";
static const std::string MKINDEX_SEGLEN_DEFAULT = "8192";
static const std::string MKINDEX_SEGLEN_DESCR   = "\
\tNumber of letters in sequence store segment. It is recommended that this \
value is set to less than half of length of a typical sequence in the \
reference database. Each reference sequence occupies at least one segment \
and the sequence store can store at most 2^32 - 1 bases. If the reference \
has a large number of very short sequence, decreasing this value can help \
to pack more sequences into the sequence store and optimize memory usage. \
The value must be a power of 2 between 32 and 8192.\n\
";

static const std::string MKINDEX_ALEXT_KEY      = "al-extend";
static const std::string MKINDEX_ALEXT_SKEY     = "";
static const std::string MKINDEX_ALEXT_LABEL    = "bases";
static const std::string MKINDEX_ALEXT_DEFAULT  = "2000";
static const std::string MKINDEX_ALEXT_DESCR    = "\
\tNumber of reference bases by which an alternate locus is extended to the \
left (right) in the case of non-fuzzy left (right) end.\n\
";

//------------------------------------------------------------------------------
static common::CFileBase::TCompression Str2Compr( const std::string & name )
{
    if(      name == "none" )  return common::CFileBase::COMPRESSION_NONE;
    else if( name == "gzip" )  return common::CFileBase::COMPRESSION_ZIP;
    else if( name == "bzip2" ) return common::CFileBase::COMPRESSION_BZIP;
    else return common::CFileBase::COMPRESSION_AUTO;
}

//------------------------------------------------------------------------------
struct CSrPrismException : public CException
{
    virtual const std::string ErrorMessage( TErrorCode ) const
    { return "srprism exception"; }

    M_EXCEPT_CTOR( CSrPrismException )
};

//------------------------------------------------------------------------------
void SetArgsForSearch( COptionsParser & options_parser ) {
    options_parser.NewGroup( "SEARCH PARAMETERS:" );
    options_parser.AddParam(
            SEARCH_INDEX_KEY, SEARCH_INDEX_SKEY,
            SEARCH_INDEX_DESCR, SEARCH_INDEX_LABEL );
    options_parser.AddParam(
            SEARCH_INPUT_KEY, SEARCH_INPUT_SKEY,
            SEARCH_INPUT_DESCR, SEARCH_INPUT_LABEL );
    options_parser.AddOptionalParam(
            SEARCH_OUTPUT_KEY, SEARCH_OUTPUT_SKEY,
            SEARCH_OUTPUT_DESCR, SEARCH_OUTPUT_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_MODE_KEY, SEARCH_MODE_SKEY, SEARCH_MODE_DEFAULT,
            SEARCH_MODE_DESCR, SEARCH_MODE_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_NERR_KEY, SEARCH_NERR_SKEY, SEARCH_NERR_DEFAULT,
            SEARCH_NERR_DESCR, SEARCH_NERR_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_DIST_KEY, SEARCH_DIST_SKEY, SEARCH_DIST_DEFAULT,
            SEARCH_DIST_DESCR, SEARCH_DIST_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_FUZZ_KEY, SEARCH_FUZZ_SKEY, SEARCH_FUZZ_DEFAULT,
            SEARCH_FUZZ_DESCR, SEARCH_FUZZ_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_INFMT_KEY, SEARCH_INFMT_SKEY, SEARCH_INFMT_DEFAULT,
            SEARCH_INFMT_DESCR, SEARCH_INFMT_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_ICOMPR_KEY, SEARCH_ICOMPR_SKEY, SEARCH_ICOMPR_DEFAULT,
            SEARCH_ICOMPR_DESCR, SEARCH_ICOMPR_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_MEM_KEY, SEARCH_MEM_SKEY, SEARCH_MEM_DEFAULT,
            SEARCH_MEM_DESCR, SEARCH_MEM_LABEL );
    options_parser.AddOptionalParam(
            SEARCH_BATCH_KEY, SEARCH_BATCH_SKEY,
            SEARCH_BATCH_DESCR, SEARCH_BATCH_LABEL );
    options_parser.AddOptionalParam(
            SEARCH_SBATCH_KEY, SEARCH_SBATCH_SKEY,
            SEARCH_SBATCH_DESCR, SEARCH_SBATCH_LABEL );
    options_parser.AddOptionalParam(
            SEARCH_EBATCH_KEY, SEARCH_EBATCH_SKEY,
            SEARCH_EBATCH_DESCR, SEARCH_EBATCH_LABEL );
    options_parser.AddFlag( 
            SEARCH_QID_KEY, SEARCH_QID_SKEY, SEARCH_QID_DESCR );
    options_parser.AddFlag( 
            SEARCH_SID_KEY, SEARCH_SID_SKEY, SEARCH_SID_DESCR );
    options_parser.AddDefaultParam(
            SEARCH_TMPDIR_KEY, SEARCH_TMPDIR_SKEY, SEARCH_TMPDIR_DEFAULT,
            SEARCH_TMPDIR_DESCR, SEARCH_TMPDIR_LABEL );
    options_parser.AddOptionalParam(
            SEARCH_PAIRED_LOG_KEY, SEARCH_PAIRED_LOG_SKEY,
            SEARCH_PAIRED_LOG_DESCR, SEARCH_PAIRED_LOG_LABEL );
    options_parser.AddParam(
            SEARCH_PAIRED_KEY, SEARCH_PAIRED_SKEY,
            SEARCH_PAIRED_DESCR, SEARCH_PAIRED_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_SKIP_UNMAPPED_KEY, SEARCH_SKIP_UNMAPPED_SKEY,
            SEARCH_SKIP_UNMAPPED_DEFAULT, SEARCH_SKIP_UNMAPPED_DESCR,
            SEARCH_SKIP_UNMAPPED_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_REPEAT_KEY, SEARCH_REPEAT_SKEY,
            SEARCH_REPEAT_DEFAULT, SEARCH_REPEAT_DESCR,
            SEARCH_REPEAT_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_RESCONF_KEY, SEARCH_RESCONF_SKEY,
            SEARCH_RESCONF_DEFAULT, SEARCH_RESCONF_DESCR,
            SEARCH_RESCONF_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_SA_START_KEY, SEARCH_SA_START_SKEY,
            SEARCH_SA_START_DEFAULT, SEARCH_SA_START_DESCR,
            SEARCH_SA_START_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_SA_END_KEY, SEARCH_SA_END_SKEY,
            SEARCH_SA_END_DEFAULT, SEARCH_SA_END_DESCR,
            SEARCH_SA_END_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_EXTRA_TAGS_KEY, SEARCH_EXTRA_TAGS_SKEY,
            SEARCH_EXTRA_TAGS_DEFAULT, SEARCH_EXTRA_TAGS_DESCR,
            SEARCH_EXTRA_TAGS_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_NRES_KEY, SEARCH_NRES_SKEY, SEARCH_NRES_DEFAULT,
            SEARCH_NRES_DESCR, SEARCH_NRES_LABEL );
    /*
    options_parser.AddDefaultParam(
            SEARCH_OUTFMT_KEY, SEARCH_OUTFMT_SKEY, SEARCH_OUTFMT_DEFAULT,
            SEARCH_OUTFMT_DESCR, SEARCH_OUTFMT_LABEL );
    */
    options_parser.AddFlag(
            SEARCH_SD_KEY, SEARCH_SD_SKEY, SEARCH_SD_DESCR );
    options_parser.AddFlag(
            SEARCH_SD_STOP_KEY, SEARCH_SD_STOP_SKEY, 
            SEARCH_SD_STOP_DESCR );
    options_parser.AddDefaultParam(
            SEARCH_SD_HNAME_KEY, SEARCH_SD_HNAME_SKEY,
            SEARCH_SD_HNAME_DEFAULT, SEARCH_SD_HNAME_DESCR,
            SEARCH_SD_HNAME_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_RANDOMIZE_KEY, SEARCH_RANDOMIZE_SKEY,
            SEARCH_RANDOMIZE_DEFAULT, SEARCH_RANDOMIZE_DESCR,
            SEARCH_RANDOMIZE_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_RANDOM_SEED_KEY, SEARCH_RANDOM_SEED_SKEY,
            SEARCH_RANDOM_SEED_DEFAULT, SEARCH_RANDOM_SEED_DESCR,
            SEARCH_RANDOM_SEED_LABEL );
    options_parser.AddDefaultParam(
            SEARCH_SAM_HEADER_KEY, SEARCH_SAM_HEADER_SKEY,
            SEARCH_SAM_HEADER_DEFAULT, SEARCH_SAM_HEADER_DESCR,
            SEARCH_SAM_HEADER_LABEL );

#ifndef NDEBUG
    options_parser.AddOptionalParam(
            SEARCH_FIX_HC_KEY, SEARCH_FIX_HC_SKEY,
            SEARCH_FIX_HC_DESCR, SEARCH_FIX_HC_LABEL );
#endif
}

//------------------------------------------------------------------------------
void SetArgsForMkIndex( COptionsParser & options_parser ) {
    options_parser.NewGroup( "MKINDEX PARAMETERS:" );
    options_parser.AddOptionalParam(
            MKINDEX_INPUT_KEY, MKINDEX_INPUT_SKEY,
            MKINDEX_INPUT_DESCR, MKINDEX_INPUT_LABEL );
    options_parser.AddOptionalParam(
            MKINDEX_INPUT_LIST_KEY, MKINDEX_INPUT_LIST_SKEY,
            MKINDEX_INPUT_LIST_DESCR, MKINDEX_INPUT_LIST_LABEL );
    options_parser.AddOptionalParam(
            MKINDEX_ALSPEC_KEY, MKINDEX_ALSPEC_SKEY,
            MKINDEX_ALSPEC_DESCR, MKINDEX_ALSPEC_LABEL );
    options_parser.AddParam(
            MKINDEX_OUTPUT_KEY, MKINDEX_OUTPUT_SKEY,
            MKINDEX_OUTPUT_DESCR, MKINDEX_OUTPUT_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_INFMT_KEY, MKINDEX_INFMT_SKEY, MKINDEX_INFMT_DEFAULT,
            MKINDEX_INFMT_DESCR, MKINDEX_INFMT_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_ICOMPR_KEY, MKINDEX_ICOMPR_SKEY, MKINDEX_ICOMPR_DEFAULT,
            MKINDEX_ICOMPR_DESCR, MKINDEX_ICOMPR_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_OUTFMT_KEY, MKINDEX_OUTFMT_SKEY, MKINDEX_OUTFMT_DEFAULT,
            MKINDEX_OUTFMT_DESCR, MKINDEX_OUTFMT_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_MEM_KEY, MKINDEX_MEM_SKEY, MKINDEX_MEM_DEFAULT,
            MKINDEX_MEM_DESCR, MKINDEX_MEM_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_SEGLEN_KEY, MKINDEX_SEGLEN_SKEY,
            MKINDEX_SEGLEN_DEFAULT, MKINDEX_SEGLEN_DESCR,
            MKINDEX_SEGLEN_LABEL );
    options_parser.AddDefaultParam(
            MKINDEX_ALEXT_KEY, MKINDEX_ALEXT_SKEY,
            MKINDEX_ALEXT_DEFAULT, MKINDEX_ALEXT_DESCR,
            MKINDEX_ALEXT_LABEL );
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // recreate command line for SAM output
    //
    static std::string CMDLINE;

    if( argc > 0 )
    {
        CMDLINE = argv[0];

        for( int i( 1 ); i < argc; ++i )
        {
            CMDLINE += " ";
            CMDLINE += argv[i];
        }
    }

    // return codes
    static const int SUCCESS   = 0;
    static const int UNKNOWN   = 1;
    static const int STDEXCEPT = 2;
    static const int EXCEPT    = 3;
    static const int NO_ARGS   = 4;

    // command strings
    static const char * SEARCH_CMD = "search";
    static const char * MKIDX_CMD  = "mkindex";
    static const char * HELP_CMD   = "help";

    static const char * HELP_PROMPT = "\n"
        "Please type 'srprism help' for general usage information;\n"
        "       type 'srprism help search' for help on 'search' command;\n"
        "       type 'srprism help mkindex' for help on 'mkindex' command.";

    std::string usage_string;

    try{
        if( argc < 2 )
        {
            std::cerr << HELP_PROMPT << std::endl;
            return NO_ARGS;
        }

        // get the srprizm command
        //
        COptionsParser options_parser( 
                argv[0], PROGRAM_DESCRIPTION, GetVersionString() );
        options_parser.NewGroup( "COMMON PARAMETERS:" );
        options_parser.AddPositionalDescription( CMD_LABEL, CMD_DESCR );
        options_parser.SetMaxOptPositionals( 1 );

        // definining common options
        //
        options_parser.AddDefaultParam(
                TRACE_LEVEL_KEY, TRACE_LEVEL_SKEY, TRACE_LEVEL_DEFAULT,
                TRACE_LEVEL_DESCR, TRACE_LEVEL_LABEL );
        options_parser.AddOptionalParam(
                LOG_FNAME_KEY, LOG_FNAME_SKEY,
                LOG_FNAME_DESCR, LOG_FNAME_LABEL );

        options_parser.Parse( argv + 1, argv + argc, false );

        std::string command;
        options_parser.Bind( CMD_POS, command );

        // setting up tracing
        //
        {
            std::string lvl_str;
            options_parser.Bind( TRACE_LEVEL_KEY, lvl_str );
            static const char * lvl_strings[] = {
                "debug", "info", "warning", "error", "quiet" };
            const char ** idx_end = lvl_strings + 5;
            const char ** idx = std::find( lvl_strings, idx_end, lvl_str );

            if( idx == idx_end ) {
                std::string allowed_strings = rfold(
                        CPlus< std::string, const char * >( " " ), 
                        std::string( "[ " ), lvl_strings, idx_end ) + "]";
                M_THROW( CSrPrismException, NOCODE,
                         "invalid trace level value " << lvl_str <<
                         "; allowed values: " << allowed_strings );
            }

            CTracer::SetLevel( idx - lvl_strings );

            if( options_parser.IsPresent( LOG_FNAME_KEY ) ) {
                std::string log_name;
                options_parser.Bind( LOG_FNAME_KEY, log_name );
                CTracer::SetOutputFile( log_name );
            }
            else CTracer::SetOutputStream( std::cerr );
        }

        if( command == SEARCH_CMD ) {
            SetArgsForSearch( options_parser );
        }
        else if( command == MKIDX_CMD ) {
            SetArgsForMkIndex( options_parser );
        }
        else if( command == HELP_CMD ) {
            std::string help_sec;

            if( HLP_SEC_POS < options_parser.NPositionals() ) {
                options_parser.Bind( HLP_SEC_POS, help_sec );
            }

            if( help_sec == SEARCH_CMD ) {
                SetArgsForSearch( options_parser );
            }
            else if( help_sec == MKIDX_CMD ) {
                SetArgsForMkIndex( options_parser );
            }

            std::cout << options_parser.Usage()
                      << std::endl;
            return 0;
        }
        else {
            M_THROW( CSrPrismException, NOCODE, 
                     "bad srprism command: " << command );
        }

        options_parser.Parse( argv + 1, argv + argc );
        seq::InitCoding();
        
        if( command == SEARCH_CMD ) {
            bool no_qids, no_sids;
            CSearch::SOptions options;
            options.cmdline = CMDLINE;
            options.force_paired = options.force_unpaired = false;
            options.start_batch = 1;
            options.end_batch = common::SIntTraits< Uint4 >::MAX;
            options.batch_limit = 10000000UL;
            options.strict_batch = false;
            options_parser.Bind( SEARCH_INDEX_KEY , options.index_basename );
            options_parser.Bind( SEARCH_NERR_KEY  , options.n_err );
            options_parser.Bind( SEARCH_DIST_KEY  , options.pair_distance );
            options_parser.Bind( SEARCH_FUZZ_KEY  , options.pair_fuzz );
            options_parser.Bind( SEARCH_INPUT_KEY , options.input );
            options_parser.Bind( SEARCH_INFMT_KEY , options.input_fmt );
            options_parser.Bind( SEARCH_MEM_KEY   , options.mem_limit );
            options_parser.Bind( SEARCH_QID_KEY   , no_qids );
            options_parser.Bind( SEARCH_SID_KEY   , no_sids );
            options_parser.Bind( SEARCH_SD_KEY    , options.discover_sep );
            options_parser.Bind( 
                    SEARCH_SD_STOP_KEY , options.discover_sep_stop );
            options_parser.Bind( SEARCH_SD_HNAME_KEY , options.hist_fname );
            options_parser.Bind( SEARCH_RANDOMIZE_KEY, options.randomize );
            options_parser.Bind( SEARCH_RANDOM_SEED_KEY, options.random_seed );
            options_parser.Bind( SEARCH_TMPDIR_KEY, options.tmpdir );
            options_parser.Bind(
                    SEARCH_REPEAT_KEY, options.repeat_threshold );
            options_parser.Bind( SEARCH_RESCONF_KEY, options.resconf_str );
            options_parser.Bind( SEARCH_SA_START_KEY, options.sa_start );
            options_parser.Bind( SEARCH_SA_END_KEY, options.sa_end );
            options_parser.Bind( SEARCH_EXTRA_TAGS_KEY, options.extra_tags );
            options_parser.Bind( SEARCH_SAM_HEADER_KEY, options.sam_header );

            {
                std::string search_mode_str;
                options_parser.Bind( SEARCH_MODE_KEY, search_mode_str );

                if( search_mode_str == "min-err" ) {
                    options.search_mode = SSearchMode::DEFAULT;
                }
                else if( search_mode_str == "bound-err" ) {
                    options.search_mode = SSearchMode::BOUND_ERR;
                }
                else if( search_mode_str == "sum-err" ) {
                    options.search_mode = SSearchMode::SUM_ERR;
                }
                else if( search_mode_str == "partial" ) {
                    options.search_mode = SSearchMode::PARTIAL;
                }
                else options.search_mode = -1;
            }

            {
                std::string compr_str;
                options_parser.Bind( SEARCH_ICOMPR_KEY, compr_str );
                options.input_compression = Str2Compr( compr_str );
            }

            if( command == SEARCH_CMD ) {
                options_parser.Bind( SEARCH_NRES_KEY  , options.res_limit );
                // options_parser.Bind( SEARCH_OUTFMT_KEY, options.output_fmt );
                options.output_fmt = "sam";
                options_parser.Bind( 
                        SEARCH_SKIP_UNMAPPED_KEY, options.skip_unmapped );
            }
            else {
                options.res_limit = 2;
                options.output_fmt = "sam";
                options.skip_unmapped = false;
            }

            options.use_qids = !no_qids;
            options.use_sids = !no_sids;

            if( options_parser.IsPresent( SEARCH_OUTPUT_KEY ) ) {
                options_parser.Bind( SEARCH_OUTPUT_KEY, options.output );
            }

            if( options_parser.IsPresent( SEARCH_PAIRED_LOG_KEY ) ) {
                std::string val;
                options_parser.Bind( SEARCH_PAIRED_LOG_KEY, val );
                options.paired_log = val;
            }
            else options.paired_log = "";

            if( options_parser.IsPresent( SEARCH_PAIRED_KEY ) ) {
                bool val;
                options_parser.Bind( SEARCH_PAIRED_KEY, val );

                if( val ) options.force_paired   = true;
                else      options.force_unpaired = true;
            }

            if( options_parser.IsPresent( SEARCH_BATCH_KEY ) ) {
                options_parser.Bind( SEARCH_BATCH_KEY, options.batch_limit );
                options.strict_batch = true;
            }

            if( options_parser.IsPresent( SEARCH_SBATCH_KEY ) ) {
                options_parser.Bind( SEARCH_SBATCH_KEY, options.start_batch );
                options.strict_batch = true;
            }

            if( options_parser.IsPresent( SEARCH_EBATCH_KEY ) ) {
                options_parser.Bind( SEARCH_EBATCH_KEY, options.end_batch );
                options.strict_batch = true;
            }

#ifndef NDEBUG
            if( options_parser.IsPresent( SEARCH_FIX_HC_KEY ) ) {
                options.use_fixed_hc = true;
                options_parser.Bind( SEARCH_FIX_HC_KEY, options.fixed_hc );
            }
#endif

            CSearch search( options );
            search.Run();
        }
        else if( command == MKIDX_CMD ) {
            CMkIdx::SOptions options;
            options_parser.Bind( MKINDEX_INFMT_KEY,  options.infmt );
            options_parser.Bind( MKINDEX_OUTFMT_KEY, options.outfmt );
            options_parser.Bind( MKINDEX_OUTPUT_KEY, options.output );
            options_parser.Bind( MKINDEX_MEM_KEY,    options.max_mem );
            options_parser.Bind( MKINDEX_SEGLEN_KEY, options.ss_seg_len );
            options_parser.Bind( MKINDEX_ALEXT_KEY,  options.al_extend );

            {
                std::string compr_str;
                options_parser.Bind( MKINDEX_ICOMPR_KEY, compr_str );
                options.input_compression = Str2Compr( compr_str );
            }

            if( options_parser.IsPresent( MKINDEX_INPUT_KEY ) ) {
                options_parser.Bind( MKINDEX_INPUT_KEY, options.input );
            }
            else if( options_parser.IsPresent( MKINDEX_INPUT_LIST_KEY ) ) {
                options_parser.Bind( 
                        MKINDEX_INPUT_LIST_KEY, options.input_list );
            }

            if( options_parser.IsPresent( MKINDEX_ALSPEC_KEY ) ) {
                options_parser.Bind(
                        MKINDEX_ALSPEC_KEY, options.alt_loc_spec );
            }

            CMkIdx mkidx( options );
            mkidx.Run();
        }
        else SRPRISM_ASSERT( false );
    }
    catch( const CException & e ) {
        M_TRACE( CTracer::ERROR_LVL, e.what() );
        std::cerr << e.what() << std::endl;
        std::cerr << HELP_PROMPT << std::endl;
        return EXCEPT;
    }
    catch( const std::exception & e ) {
        M_TRACE( CTracer::ERROR_LVL, e.what() );
        std::cerr << e.what() << std::endl;
        std::cerr << HELP_PROMPT << std::endl;
        return STDEXCEPT;
    }
    catch( ... ) {
        static const char * msg = "unknown exception";
        M_TRACE( CTracer::ERROR_LVL, msg );
        std::cerr << msg << std::endl;
        std::cerr << HELP_PROMPT << std::endl;
        return UNKNOWN;
    }

    return SUCCESS;
}

