/*  $Id: index_iterator.hpp 358393 2012-04-02 16:20:45Z morgulis $
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
 * File Description: high level code for .idx files traversal
 *
 */

#ifndef __SRPRISM_INDEX_ITERATOR_HPP__
#define __SRPRISM_INDEX_ITERATOR_HPP__

#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <string>

#include "../seq/seqdef.hpp"
#include "srprismdef.hpp"
#include "idxmap_reader.hpp"
#include "idx_reader.hpp"
#include "index_base.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <string>

#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/idxmap_reader.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/idx_reader.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/index_base.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CIndexIterator : public CIndexBase
{
    static const size_t LETTER_BITS = 
        seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;

    static const unsigned int COMPRESSION = common::BYTEBITS/LETTER_BITS;

    static const size_t PREFIX_BYTES = 1 + (MAP_PREFIX_LEN - 1)/COMPRESSION;
    static const size_t SFX_BITS = LETTER_BITS*(PREFIX_LEN - MAP_PREFIX_LEN);
    static const size_t EXACT_BYTES = 1 + (3 + SFX_BITS)/common::BYTEBITS;
    
    public:

        typedef common::Uint4 TPos;
        typedef common::Uint8 TUnit;

    private:

        static const size_t EXT_BITS = EXT_LEN*LETTER_BITS;
        static const size_t EXT_SHIFT = 
            common::BYTEBITS*sizeof( TUnit ) - EXT_BITS;

        typedef std::vector< TPos > TPosVec;

        struct SExtInfo
        {
            TWord Extension( void ) const { return (TWord)ext; }

            TUnit ext;
            size_t fnpos, rnpos, index;
        };

    public:

        typedef SExtInfo TExtInfo;
        typedef std::vector< SExtInfo > TExtData;
        typedef TPosVec::const_iterator TPosIter;

        CIndexIterator( const std::string & basename );

        bool Seek( TUnit prefix );
        TUnit Prefix(void) const { return prefix_; }
        bool Special(void) const { return special_; }

        size_t NPos( void ) const
        { return Special() ? pos_.size()/2 : pos_.size(); }

        TPosIter PosStart( TStrand strand ) const 
        { 
            return (strand == seq::STRAND_FW || palindrome_ ) 
                   ? pos_.begin() 
                   : pos_.begin() + rv_pos_start_; 
        }

        TPosIter PosEnd( TStrand strand ) const
        {
            return (strand == seq::STRAND_RV || palindrome_ ) 
                   ? pos_.end()
                   : pos_.begin() + rv_pos_start_;
        }

        const TExtData & Extensions( TStrand strand ) const
        { return ext_data_[strand]; }

    private:

        CIndexIterator( const CIndexIterator & );
        CIndexIterator & operator=( const CIndexIterator & );

        static CIdxMapReader::TUnit Unit2Map( TUnit unit )
        { return (CIdxMapReader::TUnit)(unit>>SFX_BITS); }

        static TUnit Map2Unit( CIdxMapReader::TUnit v )
        { return (v<<SFX_BITS); }

        bool Next(void);
        void ReadDataSpecial(void);
        void ReadData(void);
        template< TStrand strand > void ReadStrandDataSpecial(void);
        void CheckPalindrome(void);

        CIdxMapReader map_reader_;
        CIdxReader idx_reader_;

        enum EState
        {
            NEW_PREFIX,
            IN_PREFIX,
            END_OF_INDEX
        } state_;

        TUnit prefix_;
        size_t bytes_left_;
        bool init_;
        bool special_;
        bool palindrome_;
        bool end_;
        bool start_;
        size_t r_bytes_, f_bytes_;
        size_t rv_pos_start_;
        TExtData ext_data_[seq::N_STRANDS];
        TPosVec pos_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

