/*  $Id: srprismdef.hpp 590234 2019-07-25 16:28:25Z morgulis $
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
 * Authors:  Aleksandr Morgulis
 *
 * File Description: common definitions for srprism 
 *
 */

#ifndef __SRPRISM__SRPRISMDEF_HPP__
#define __SRPRISM__SRPRISMDEF_HPP__

#include "../common/def.h"

#include <sstream>

#ifndef NCBI_CPP_TK

#include <common/bits.hpp>
#include <seq/seqdef.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
// ALGORITHM VERSION INFORMATION
//
static const int VERSION_MAJOR = 3;
static const int VERSION_MINOR = 1;
static const int VERSION_POINT = 0;
static const int VERSION_STATE = 0;

inline std::string GetVersionString( void )
{
    std::ostringstream os;
    os << VERSION_MAJOR << '.' << VERSION_MINOR << '.' << VERSION_POINT;

    switch( VERSION_STATE ) {
        case 1: os << "-beta"; break;
        case 2: os << "-alpha"; break;
    }
    
    return os.str();
}

//------------------------------------------------------------------------------
//  Some types used throughout srprism.
//
typedef common::Uint8 TSize;

//------------------------------------------------------------------------------
/** 
  \defgroup forward Forwarding Declarations
 */
/*@{*/
typedef seq::TSeqSize TSeqSize;
typedef seq::TLetter  TLetter;
typedef seq::TStrand  TStrand;
typedef seq::TCoding  TCoding;
/*@}*/

//------------------------------------------------------------------------------
/** Maximum supported sequence length. */
static const TSeqSize MAX_SEQ_SIZE = common::SIntTraits< TSeqSize >::MAX;

/** Representation of an invelid sequence position. */
static const TSeqSize INVALID_SEQ_POS = MAX_SEQ_SIZE;

//------------------------------------------------------------------------------
/**
  \defgroup aligndef Declarations Related to Alignments
  */
/*@{*/
typedef common::Uint1 TErrType; /**< Possible error combinations. */

/** Encoding of letters in alignment errors. */
static const TCoding ERROR_CODING = seq::CODING_NCBI4NA;

/** Letter representing a gap in the alignement. */
static const TLetter ERROR_GAP_LETTER = 
    seq::SCodingTraits< ERROR_CODING >::GAP_LETTER;

//------------------------------------------------------------------------------
/** Symbolic names for different error combinations. */
struct SErrType {
    static const TErrType E = 0; /**< Exact match. */
    static const TErrType M = 1; /**< Mismatch. */
    static const TErrType I = 2; /**< Insertion in the query. */
    static const TErrType D = 3; /**< Deletion from the query. */

    /** Maximum error type value describing single errors. */
    static const TErrType ONE_ERROR_MAX_VAL = D;

    static const TErrType MM = 4;  /**< Double mismatch. */
    static const TErrType MI = 5;  /**< Mismatch followed by insertion. */
    static const TErrType MD = 6;  /**< Mismatch followed by deletion. */
    static const TErrType IM = 7;  /**< Insertion followed by mismatch. */
    static const TErrType DM = 8;  /**< Deletion followed by mismatch. */
    static const TErrType II = 9;  /**< Double insertion. */
    static const TErrType DD = 10; /**< Double deletion. */
    static const TErrType ID = 11; /**< Insertion followed by deletion. */
    static const TErrType DI = 12; /**< Deletion followed by insertion. */

    /** Maximum error type value describing double errors. */
    static const TErrType TWO_ERROR_MAX_VAL = DI;

    /** Total number of different error combination types. */
    static const TErrType NUM_ERR_TYPES = TWO_ERROR_MAX_VAL + 1;

    /** Get the number of errors from the given error type.

        \param t error type
        \return number of errors in t
    */
    static int NErr( TErrType t )
    { 
        SRPRISM_ASSERT( t <= TWO_ERROR_MAX_VAL );
        return (t > ONE_ERROR_MAX_VAL) ? 2 : ((t > E) ? 1 : 0); 
    }

    /** Get error type based on alignment letters.
      
        Note: does not work with CODING == CODING_NCBI2NA.

        \param ql query letter
        \param sl subject letter

        \return error type
      */
    template< TCoding CODING >
    static TErrType ErrType( TLetter ql, TLetter sl )
    {
        static const TLetter GAP_LETTER = 
            seq::SCodingTraits< CODING >::GAP_LETTER; 

        if( ql == GAP_LETTER ) {
            if( sl == GAP_LETTER ) return E;
            else return D;
        }
        else if( sl == GAP_LETTER ) return I;
        else return M;
    }
};
/*@}*/

//------------------------------------------------------------------------------
/** Max number of errors supported by srprism search. */
const int MAX_ERR = 15;

/** Max number of errors for short queries. */
const common::Uint1 MAX_SHORT_ERR = 1;

/** Max number of errors in the medium query seed. */
const common::Uint1 MAX_MED_SEED_N_ERR = 2;

/** Max number of errors in a seeding area */
const common::Uint1 MAX_SEED_N_ERR = 15;

/** Max number of hashes in a query seeding area. */
const int MAX_N_HASHES = MAX_SEED_N_ERR + 1;

/** Main encoding used for sequence data. */
const TCoding SEQDATA_CODING = seq::CODING_NCBI2NA;

/** Encoding used to output sequence data. */
const TCoding OUTPUT_CODING = seq::CODING_IUPACNA;

/** Type to represent query sequence position within a query stream. */
typedef common::Uint8 TQueryOrdId;

/** Type to represent subject sequence position within a database. */
typedef common::Uint4 TDBOrdId;

/** Type to represent query sequence position within a batch. */
typedef common::Uint4 TQNum;

/** Size type for storing number of results per query. */
typedef common::Uint4 TNumRes;

/** Nmer size used for query hashing */
const TSeqSize HASH_LEN = 16;

/** Length of the prefix on which index is sorted. */
const TSeqSize PREFIX_LEN = HASH_LEN;

/** Length of the extension for bad n-mers processing. */
const TSeqSize EXT_LEN = HASH_LEN;

/** Lower limit on query length. */
const TSeqSize MIN_QUERY_LEN = PREFIX_LEN;

/** Lower limit on `long` queries.

    Queries shorter than this are searched with at most 1 error.
*/
const TSeqSize MIN_MED_QUERY_LEN = PREFIX_LEN + EXT_LEN;

/** Max length of the query seed area */
const TSeqSize MIN_LONG_QUERY_LEN = 3*HASH_LEN;

/** Upper imit on query length. */
const TSeqSize MAX_QUERY_LEN = (TSeqSize)(8*common::KILOBYTE);

/** Integer type used to represent 16-mer prefixes. */
typedef common::Uint4 TPrefix;

/** Integer type used to represent 16-mer extensions. */
typedef common::Uint4 TExtension;

/** Word and half word types used for internal subject and query storage. */
typedef common::Uint4 TWord;
typedef common::Uint2 THalfWord;

//------------------------------------------------------------------------------
// FILE NAMES
//
extern const char * QDUMP_NAME; /**< temporary containing batch query data */

//------------------------------------------------------------------------------
// PAIRED RESULT CONFIGURATIONS
//
// IPAM = in-place alignment mask
//
// IPAM value is a 4-bit integer that specifies how to search for a
// mate alignment when the initial alignment is obtained. Four possible
// configurations are encoded as:
//
// 0 - search left interval; forward seeds;
// 1 - search left interval; reverse seeds;
// 2 - search right interval; forward seeds;
// 3 - search right interval; reverse seeds;
//
// In IPAM value a bit with index i is on if the i-th configuration should
// be searched.
//
// An IPAM vector contains 4 IPAM values corresponding to the following
// 4 cases:
//
// 0 - initial alignment is for the first mate and is forward-aligned;
// 1 - initial alignment is for the first mate and is reverse-aligned;
// 2 - initial alignment is for the second mate and is forward-aligned;
// 3 - initial alignment is for the second mate and is reverse-aligned;
//
static const size_t MAX_IPAM_IDX = 3; // size of IPAM vector
typedef common::Uint1 T_IPAM;         // type for IPAM values
struct S_IPAM { T_IPAM data[MAX_IPAM_IDX + 1]; }; // type of IPAM vector

// masks to determine which subject intervals to search for mate alignment
static const T_IPAM IPAM_LEFT_ENABLED  = (T_IPAM)0x3;
static const T_IPAM IPAM_RIGHT_ENABLED = (T_IPAM)0xc;

// masks to determine which direction to search for mate alignment
static const T_IPAM IPAM_FW_ENABLED = (T_IPAM)0x5;
static const T_IPAM IPAM_RV_ENABLED = (T_IPAM)0xa;

//------------------------------------------------------------------------------
// NAMES FOR DIFFERENT PASS TYPES
//
static const int HASH_NORMAL = 0; // pass processing hashes 0, 1, and 2
static const int HASH_BLOWUP = 1; // pass processing hash 3 (blowup)

static const size_t N_HASHES = 3; // max number of normal hashes in a seeding
                                  //    area

extern S_IPAM ParseResConfStr( const std::string rcstr );

//------------------------------------------------------------------------------
// COUNTER NAMES FOR SEARCH STATISTICS
//
extern const char * STAT_N_ALIGNS;
extern const char * STAT_N_UALIGNS;
extern const char * STAT_N_FILTER;
extern const char * STAT_N_CANDIDATES;
extern const char * STAT_N_INPLACE;
extern const char * STAT_N_INPLACE_ALIGNS;

END_NS( srprism )
END_STD_SCOPES

#endif

