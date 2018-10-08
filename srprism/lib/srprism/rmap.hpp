/*  $Id: rmap.hpp 219989 2011-01-14 21:41:34Z morgulis $
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
 * File Description: repeat map interface
 *
 */

#ifndef __SRPRISM_RMAP_HPP__
#define __SRPRISM_RMAP_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <vector>
#include <string>
#include <algorithm>

#include <srprism/srprismdef.hpp>
#include <srprism/index_base.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <vector>
#include <string>
#include <algorithm>

#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/index_base.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
//
// the structure that maps 16-mers to their repeat counts in the genome;
// the size is kept in check by not keeping data for 16-mers with 
// binlog( repeat count ) < NPOS_LOG_START;
// binary logs of the counts are stores, rather than counts themselves;
// lookup is logarithmic via a binary search on 16-mer values
//
class CRMap : public CIndexBase
{
    private:

        // 16-mer list type
        //
        typedef std::vector< TPrefix > TPrefixes;

        // counts list type
        //
        typedef std::vector< common::Uint1 > TRanks;

    public:

        // instance constructor from the index base name
        //
        CRMap( const std::string & name );

        // get repeat count (logarithm) by 16-mer value
        //
        common::Uint1 RepeatRank( TPrefix prefix ) const
        {
            TPrefixes::const_iterator i( std::lower_bound( 
                        prefixes_.begin(), prefixes_.end(), prefix ) );

            if( i == prefixes_.end() || *i != prefix ) return 0;
            else return ranks_[i - prefixes_.begin()];
        }

    private:

        TPrefixes prefixes_;    // list of kept 16-mers
        TRanks ranks_;          // list of corresponding counts (logarithms)
};

END_NS( srprism )
END_STD_SCOPES

#endif

