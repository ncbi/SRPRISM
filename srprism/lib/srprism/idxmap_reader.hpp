/*  $Id: idxmap_reader.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: low-level reading and buffering of .map files
 *
 */

#ifndef __SRPRISM_IDXMAP_READER_HPP__
#define __SRPRISM_IDXMAP_READER_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <string>
#include <vector>

#include <common/exception.hpp>
#include <common/binfile.hpp>
#include <seq/seqdef.hpp>
#include <srprism/srprismdef.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <string>
#include <vector>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/binfile.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CIdxMapReader
{
    static const size_t BUFSIZE = (size_t)(4*1024);

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode SIZE = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == SIZE ) return "size error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        static const TSeqSize MAP_PREFIX_LEN = 12;

        typedef size_t TOffset;
        typedef common::Uint4 TUnit;

        CIdxMapReader( const std::string & name );

        bool Seek( TUnit target );
        bool Advance(void) { return Seek( start_ + curr_ + 1 ); }
        TOffset Offset() const { return buf_[curr_]; }
        TUnit Unit() const { return (TUnit)(start_ + curr_); }

    private:

        CIdxMapReader( const CIdxMapReader & );
        CIdxMapReader & operator=( const CIdxMapReader & );

        static const TOffset END_VAL = (TOffset)0xFFFFFFFFFFFFFFFFULL;

        static const seq::TCoding CODING = SEQDATA_CODING;
        static const size_t NUM_UNITS = 
            common::SBitFieldTraits< 
                size_t, 
                MAP_PREFIX_LEN*seq::SCodingTraits< 
                    CODING >::LETTER_BITS >::MAX + 1;

        typedef std::vector< TOffset > TBuf;

        bool NextBuf(void);

        TBuf buf_;
        common::CReadBinFile is_;

        size_t sz_;
        size_t start_;
        size_t curr_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

