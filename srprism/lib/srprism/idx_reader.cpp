/*  $Id: idx_reader.cpp 639133 2021-10-13 17:24:59Z morgulis $
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

#include <ncbi_pch.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cstring>

#ifndef WIN32
#	include <unistd.h>

#	define OPEN  open
#	define LSEEK lseek
#	define READ  read

#   define OPEN_FLAGS (O_RDONLY)
#else
#	include <io.h>

#	define OPEN  _open
#	define LSEEK _lseeki64
#	define READ  _read

#   define OPEN_FLAGS (O_RDONLY|O_BINARY)
#endif

#include <errno.h>
#include <cstring>

#include "../common/util.hpp"
#include "../common/trace.hpp"
#include "idx_reader.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
namespace {
    static const TSeqSize HEADER_LEN = 4;

    template< typename int_t >
    void ReadInt( int_t & val, int fd, size_t & off )
    {
        ssize_t bytes = READ( fd, (void *)&val, sizeof( int_t ) );

        if( bytes < 0 ) {
            M_THROW( CIdxReader::CException, SYSTEM, 
                     "read() failed [" << errno << "]: " 
                     << strerror( errno ) );
        }

        if( bytes != sizeof( int_t ) ) {
            M_THROW( CIdxReader::CException, READ,
                     "failed to read " << sizeof( int_t ) << " bytes" );
        }

        off += sizeof( int_t );
    }
}

//------------------------------------------------------------------------------
CIdxReader::CIdxReader( const std::string & name )
    : buf_( BUFSIZE, 0 ), fd_( -1 ), sz_( 0 ), start_( 0 ), off_( 0 ), 
      eof_( false ), n_seeks_( 0 ), n_reads_( 0 )
{
    fd_ = ::OPEN( name.c_str(), OPEN_FLAGS );
    
    if( fd_ < 0 ) {
        M_THROW( CException, SYSTEM, 
                 "open() failed for " << name << "; [" << errno << "]: " << 
                 strerror( errno ) );
    }

    // check endianness match and skip header
    {
        Uint1 endianness;
        ReadInt( endianness, fd_, start_ );

        if( Endianness() != (SEndianness::TTYPE)endianness ) {
            M_THROW( CException, ENDIAN, "in " << name );
        }

        for( size_t i = 1; i < HEADER_LEN; ++i ) {
            ReadInt( endianness, fd_, start_ );
        }
    }

    M_TRACE( CTracer::INFO_LVL, "opened index " << name );
}

//------------------------------------------------------------------------------
void CIdxReader::FF( TOffset offset )
{
    if( offset <= start_ + off_ ) return;

    if( offset >= start_ + sz_ ) {
        if( LSEEK( fd_, offset, SEEK_SET ) == (off_t)(-1) ) {
            M_THROW( CException, SYSTEM, 
                    "lseek() failed [" << errno << "]: " 
                    << strerror( errno ) );
        }

        ++n_seeks_;
        start_ = offset;
        sz_ = off_ = 0;
    }
    else off_ = offset - start_;
}

//------------------------------------------------------------------------------
void CIdxReader::ReadBuf(void)
{
    size_t rest = sz_ - off_;
    for( TBuf::size_type i = 0; i < rest; ++i ) buf_[i] = buf_[off_ + i];
    ssize_t n_bytes = ::READ( fd_, &buf_[0] + rest, BUFSIZE - rest );
    ++n_reads_;
    
    if( n_bytes < 0 ) {
        M_THROW( CException, SYSTEM, 
                 "read() failed [" << errno << "]: " << strerror( errno ) );
    }

    if( n_bytes == 0 ) eof_ = true;
    start_ += off_;
    off_ = 0;
    sz_ = n_bytes + rest;
}

END_NS( srprism )
END_STD_SCOPES

