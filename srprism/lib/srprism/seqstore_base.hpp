/*  $Id: seqstore_base.hpp 274326 2011-04-13 18:10:45Z morgulis $
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
 * File Description: common functionality and definitions for seqstore 
 *                   creator and reader classes
 *
 */

#ifndef __SRPRISM_SEQSTORE_BASE_HPP__
#define __SRPRISM_SEQSTORE_BASE_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <srprism/srprismdef.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CSeqStoreBase
{
    public:

        typedef common::Uint4 TPos;

    protected:

        typedef common::Uint8 TMaskUnit;
        typedef std::vector< TPos > TPosMap;
        typedef std::vector< TSeqSize > TLenMap;
        typedef std::vector< TMaskUnit > TAmbigMask;

        struct SAmbigRun
        {
            TPos pos;
            seq::TSeqSize len;
            size_t offset;

            SAmbigRun( TPos p, size_t o ) : pos( p ), len( 0 ), offset( o ) {}

            friend bool operator<( 
                    const SAmbigRun & l, const SAmbigRun & r )
            { return l.pos < r.pos; }
        };

        static const size_t SS_VERSION = 2;

        static const size_t DATA_SIZE_LIMIT = 
            1ULL + (size_t)common::SIntTraits< common::Uint4 >::MAX;
        static const size_t MASK_UNIT_BITS = 
            common::BYTEBITS*sizeof( TMaskUnit );
        static const size_t MASK_SHIFT = 
            common::SBinLog< MASK_UNIT_BITS >::VALUE;
        static const size_t MASK_MASK = 
            common::SBitFieldTraits< size_t, MASK_SHIFT >::MASK;
        static const size_t WORD_LETTERS = 
            sizeof( TWord )*seq::SCodingTraits< SEQDATA_CODING >::PACK_FACTOR;
        static const size_t WORD_BITS = sizeof( TWord )*common::BYTEBITS;
        static const size_t WORD_SHIFT = common::SBinLog< WORD_LETTERS>::VALUE;
        static const size_t WORD_MASK = 
            common::SBitFieldTraits< size_t, WORD_SHIFT >::MASK;

        static size_t MaskBits( size_t segment_letters )
        { return DATA_SIZE_LIMIT/segment_letters; }

        static size_t MaskSize( size_t segment_letters )
        { return MaskBits( segment_letters )/common::BYTEBITS; }

        static size_t MaskUnits( size_t segment_letters )
        { return MaskSize( segment_letters )/sizeof( TMaskUnit ); }

        static size_t SegmentWords( size_t segment_letters )
        { return (segment_letters>>WORD_SHIFT); }

        static size_t SegmentShift( size_t segment_letters )
        { return common::BinLog( segment_letters ); }

    public:

        static const char * DESC_SFX;
        static const char * ID_MAP_SFX;
        static const char * POS_MAP_SFX;
        static const char * SEQ_DATA_SFX;
        static const char * MASK_DATA_SFX;
        static const char * AMBIG_MAP_SFX;
        static const char * AMBIG_DATA_SFX;
};

END_NS( srprism )
END_STD_SCOPES

#endif

