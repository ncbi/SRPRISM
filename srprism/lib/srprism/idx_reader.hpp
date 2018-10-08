/*  $Id: idx_reader.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: low-level reading and buffering of .idx files
 *
 */

#ifndef __SRPRISM_IDX_READER_HPP__
#define __SRPRISM_IDX_READER_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <string>
#include <vector>
#include <algorithm>

#include <common/exception.hpp>
#include <common/trace.hpp>
#include <srprism/idxmap_reader.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <string>
#include <vector>
#include <algorithm>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/trace.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/idxmap_reader.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CIdxReader
{
    static const size_t BUFSIZE = (size_t)(32*common::KILOBYTE);

    typedef CIdxMapReader::TOffset TOffset;
    typedef std::vector< char > TBuf;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode SYSTEM = 0;
            static const TErrorCode READ   = 1;
            static const TErrorCode ENDIAN = 2;
            static const TErrorCode SIZE   = 3;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == SYSTEM ) return "system error";
                else if( code == READ ) return "read error";
                else if( code == ENDIAN ) return "endianness mismatch";
                else if( code == SIZE ) return "size error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException );
        };

        CIdxReader( const std::string & name );

        ~CIdxReader(void)
        {
            M_TRACE( common::CTracer::INFO_LVL,
                     "closing index; " << n_seeks_ << " seeks; "
                     << n_reads_ << " reads" );
        }

        void FF( TOffset offset );
        void ReadBuf(void);

        template< typename int_t >
        void NextWord( int_t & val )
        {
            if( off_ + sizeof( int_t ) > sz_ ) ReadBuf();

            if( off_ + sizeof( int_t ) > sz_ ) {
                M_THROW( CException, SIZE,
                         "unexpected end of index at position " << start_ );
            }

            TBuf::iterator i = buf_.begin() + off_;
            std::copy( i, i + sizeof( int_t ), (char *)&val );
            off_ += sizeof( int_t );
        }

        bool Eof(void) { return eof_; }

    private:

        CIdxReader( const CIdxReader & );
        CIdxReader & operator=( const CIdxReader & );

        TBuf buf_;
        int fd_;
        size_t sz_, start_, off_;
        bool eof_;
        common::Uint8 n_seeks_, n_reads_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

