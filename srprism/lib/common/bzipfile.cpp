/*  $Id: bzipfile.cpp 205312 2010-09-16 19:22:21Z morgulis $
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
 * File Description: support for bzip2 based text file compression
 *
 */

#include <errno.h>
#include <cassert>
#include <cstring>

#include <string>

#include <common/bzipfile.hpp>

START_STD_SCOPES
START_NS( common )

#ifdef WIN32
//------------------------------------------------------------------------------
CReadTextFile_BZip::CReadTextFile_BZip(const std::string& name)
    : CReadTextFile(name), eof_(false), f_(0)
{
    throw std::runtime_error("unimplemented");
}

const std::string CReadTextFile_BZip::GetLine_Impl(void)
{
    throw std::runtime_error("unimplemented");
    return "";
}

#else
//------------------------------------------------------------------------------
CReadTextFile_BZip::CReadTextFile_BZip( const std::string & name )
    : CReadTextFile( name ), eof_( false ), f_( 0 ), bzf_( 0 )
{
    f_ = fopen( name.c_str(), "r" );
    
    if( f_ == 0 ) {
        std::string errmsg( strerror( errno ) );
        M_THROW( CException, IO_ERROR, 
                 errmsg << " [" << errno << "] at opening of " << name_ );
    }

    int bzerr( 0 );
    bzf_ = BZ2_bzReadOpen( &bzerr, f_, 0, 0, NULL, 0 );

    if( bzerr != BZ_OK ) {
        if( bzerr == BZ_IO_ERROR ) {
            std::string errmsg( strerror( errno ) );
            M_THROW( CException, IO_ERROR, 
                     errmsg << " [" << errno << "] at opening of " << name_ );
        }
        else M_THROW( CException, BZIP_ERROR, 
                      "bzip error [" << bzerr << "] at opening of " << name_ );
    }
}

//------------------------------------------------------------------------------
const std::string CReadTextFile_BZip::GetLine_Impl( void )
{
    std::string result;

    while( true ) {
        char c( 0 );
        int bzerr( 0 );
        int n_read( BZ2_bzRead( &bzerr, bzf_, (void *)&c, 1 ) );

        if( bzerr == BZ_OK || bzerr == BZ_STREAM_END ) {
            if( n_read == 1 ) {
                if( c == '\n' ) { 
                    if( bzerr == BZ_STREAM_END ) eof_ = true;
                    return result; 
                }
                else if( bzerr == BZ_STREAM_END ) {
                    eof_ = true;
                    M_THROW( CException, NO_EOL, 
                             "at " << name_ << ":" << lineno_ );
                }
                else { result.append( 1, c ); }
            }
            else if( n_read == 0 ) {
                if( bzerr == BZ_STREAM_END ) {
                    eof_ = true;
                    return result;
                }
                else SRPRISM_ASSERT( false );
            }
            else SRPRISM_ASSERT( false );
        }
        else {
            if( bzerr == BZ_IO_ERROR ) {
                std::string errmsg( strerror( errno ) );
                M_THROW( CException, IO_ERROR, 
                         errmsg << " [" << errno << "] at " << 
                         name_ << ":" << lineno_ );
            }
            else M_THROW( CException, BZIP_ERROR, 
                          "bzip error [" << bzerr << "] at " << 
                          name_ << ":" << lineno_ );
        }
    }
}
#endif

END_NS( common )
END_STD_SCOPES

