/*  $Id: binfile.cpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: utilities for binary file i/o.
 *
 */

#include <ncbi_pch.hpp>

#include <stdexcept>
#include "../common/def.h"
#include "binfile.hpp"

START_STD_SCOPES
START_NS( common )

#define CHECK_STREAM(s,code,msg) if( !s.good() && !s.eof() ) {\
    M_THROW( CFileBase::CException,code,msg ); }

//------------------------------------------------------------------------------
CReadBinFile::CReadBinFile( const std::string & name )
    : CFileBase( name )
{
    try { 
        is_.open( name_.c_str(), std::ios::binary ); 
		CHECK_STREAM( is_, OPEN, "[" << name_ << "]" );
    }
    catch( std::exception & e ) {
        M_THROW( CException, SYSTEM, 
                 " at open of " << name_ <<": " << e.what() );
    }
}

//------------------------------------------------------------------------------
CReadBinFile::TSize CReadBinFile::Read( char * buf, TSize n, bool strict )
{
    try {
        is_.read( buf, n );
        CHECK_STREAM( is_, READ, "[" << name_ << ":" << pos_ << "]" );
        TSize t = is_.gcount();

        if( strict && t != n ) {
            if( t != 0 || !is_.eof() ) {
                M_THROW( CException, SIZE,
                         "failed to read record of length " << n <<
                         " at position " << pos_ );
            }
        }

        pos_ += t;
        return t;
    }
    catch( std::exception & e ) {
        M_THROW( CException, SYSTEM, 
                 " at read of " << name_ <<":" << pos_ << ": " << e.what() );
    }
}

//------------------------------------------------------------------------------
CWriteBinFile::CWriteBinFile( const std::string & name )
    : CFileBase( name )
{
    try {
        os_.open( name_.c_str(), std::ios::binary );
        if( !os_.good() ) M_THROW( CException, OPEN, "[" << name_ << "]" );
    }
    catch( std::exception & e ) {
        M_THROW( CException, SYSTEM, 
                 " at open of " << name_ << ": " << e.what() );
    }
}

//------------------------------------------------------------------------------
void CWriteBinFile::Write( const char * buf, TSize n )
{
    try {
        os_.write( buf, n );

        if( !os_.good() ) {
            M_THROW( CException, WRITE, "[" << name_ << ":" << pos_ << "]" );
        }

        pos_ += n;
    }
    catch( std::exception & e ) {
        M_THROW( CException, SYSTEM,
                 " at write of " << n << " bytes into " << name_ 
                 << ":" << pos_ << ": " << e.what() );
    }
}

#undef CHECK_STREAM

END_NS( common )
END_STD_SCOPES

