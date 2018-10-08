/*  $Id: index_base.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: types and constants used for index manipulation
 *
 */

#ifndef __SRPRISM_INDEX_BASE_HPP__
#define __SRPRISM_INDEX_BASE_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>
#include <srprism/idxmap_reader.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>
#include <../src/internal/align_toolbox/srprism/lib/srprism/idxmap_reader.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CIndexBase
{
    protected:

        static const size_t HEADER_SIZE = 4;

        static const char * IDX_MAP_SFX;
        static const char * IDX_PROPER_SFX;
        static const char * IDX_REPMAP_SFX;

        static const TSeqSize MAP_PREFIX_LEN = CIdxMapReader::MAP_PREFIX_LEN;

        static const size_t R_BYTES_START_BIT = 0;
        static const size_t R_BYTES_END_BIT   = 2;
        static const size_t F_BYTES_START_BIT = 2;
        static const size_t F_BYTES_END_BIT   = 4;
        static const size_t SFX_START_BIT     = 4;

        static const size_t NPOS_LOG_START = 7;

        struct SRMapEntry
        {
            TPrefix prefix;
            common::Uint2 flog;
            common::Uint2 rlog;
        };
};

END_NS( srprism )
END_STD_SCOPES

#endif

