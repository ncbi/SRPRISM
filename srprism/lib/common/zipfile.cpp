/*  $Id: zipfile.cpp 205312 2010-09-16 19:22:21Z morgulis $
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
 * File Description: support for zlib based text file compression
 *
 */

#include <errno.h>
#include <cstring>

#include <string>
#include <vector>

#include <common/zipfile.hpp>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
CReadTextFile_Zip::CReadTextFile_Zip( const std::string & name )
    : CReadTextFile( name )
{
    gzf_ = gzopen( name.c_str(), "r" );

    if( gzf_ == 0 ) {
        if( errno == 0 ) {
            M_THROW( CException, ZIP_ERROR, "at opening of " << name_ );
        }
        else {
            std::string errmsg( strerror( errno ) );
            M_THROW( CException, IO_ERROR, 
                     errmsg << " [" << errno << "] at opening of " << name_ );
        }
    }
}

//------------------------------------------------------------------------------
const std::string CReadTextFile_Zip::GetLine_Impl( void )
{
    std::string result;
    std::vector< char > chars( 1024, 0 );

    while( true ) {
        if( gzgets( gzf_, (char *)&chars[0], 1024 ) == Z_NULL ) {
            int errcode;
            std::string errmsg( gzerror( gzf_, &errcode ) );
            M_THROW( CException, ZIP_ERROR, 
                     "read error: " << errmsg << " [" << errcode << "]" << 
                     " at " << name_ << ":" << lineno_ );
        }

        size_t l( strlen( &chars[0] ) );
        bool eol( chars[l-1] == '\n' );

        if( eol ) {
            chars[--l] = 0;
            result.append( &chars[0] );
            return result;
        }
        else if( Eof() ) {
            M_THROW( CException, NO_EOL, "at " << name_ << ":" << lineno_ );
        }
        else result.append( &chars[0] );
    }
}

END_NS( common )
END_STD_SCOPES

