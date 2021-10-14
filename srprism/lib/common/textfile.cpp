/*  $Id: textfile.cpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: utilities for text based file i/o
 *
 */

#include <ncbi_pch.hpp>

#include "def.h"

#include <cassert>
#include <stdexcept>
#include <fstream>

#include "textfile.hpp"
#include "zipfile.hpp"
#include "bzipfile.hpp"

START_STD_SCOPES
START_NS( common )

#define CHECK_STREAM(s,code,msg) if( !s.good() && !s.eof() ) {\
    M_THROW( CFileBase::CException,code,msg ); }

//------------------------------------------------------------------------------
CReadTextFile_CPPStream::CReadTextFile_CPPStream( const std::string & name )
    : CReadTextFile( name ),
      is_( name.empty() ? std::cin : *new std::ifstream( name.c_str() ) ),
      is_holder_( name.empty() ? 0 : &is_ )
{
    try { CHECK_STREAM( is_, OPEN, "" ); }
    catch( std::exception & e ) {
        M_THROW( CFileBase::CException, SYSTEM, 
                 "at open of " << name_ << "[" << e.what() << "]" );
    }
}

//------------------------------------------------------------------------------
const std::string CReadTextFile_CPPStream::GetLine_Impl()
{
    try {
        std::string result;
        getline( is_, result );
        CHECK_STREAM( is_, READ, "" );

        if( is_.eof() && is_.gcount() > 0 ) {
            M_THROW( CException, NO_EOL, "at line " << lineno_ );
        }

        return result;
    }
    catch( CException & ) { throw; }
    catch( std::exception & e ) {
        M_THROW( CFileBase::CException, SYSTEM,
                 "at read of " << name_ << ":" << lineno_ << "(" <<
				 e.what() << ")" );
    }
}

//------------------------------------------------------------------------------
// std::auto_ptr< CReadTextFile > CReadTextFile::MakeReadTextFile( 
std::unique_ptr< CReadTextFile > CReadTextFile::MakeReadTextFile( 
        const std::string & name, TCompression c )
{
    if( c == COMPRESSION_AUTO ) c = Name2CompressionType( name );

    switch( c ) {
        case COMPRESSION_NONE:
            // return std::auto_ptr< CReadTextFile >( 
            return std::unique_ptr< CReadTextFile >( 
                    new CReadTextFile_CPPStream( name ) );
        case COMPRESSION_ZIP:
            // return std::auto_ptr< CReadTextFile >( 
            return std::unique_ptr< CReadTextFile >(
                    new CReadTextFile_Zip( name ) );
        case COMPRESSION_BZIP:
            // return std::auto_ptr< CReadTextFile >( 
            return std::unique_ptr< CReadTextFile >( 
                    new CReadTextFile_BZip( name ) );
        default: SRPRISM_ASSERT( false );
    }

    return std::unique_ptr< CReadTextFile >( nullptr );
    // return std::auto_ptr< CReadTextFile >( 0 );
}

//------------------------------------------------------------------------------
CWriteTextFile_CPPStream::CWriteTextFile_CPPStream( const std::string & name )
    : CWriteTextFile( name ),
      os_( name.empty() ? std::cout : *new std::ofstream( name.c_str() ) ),
      os_holder_( name.empty() ? 0 : &os_ )
{ if( !os_.good() ) M_THROW( CFileBase::CException, OPEN, "for " << name_ ); }

//------------------------------------------------------------------------------
// std::auto_ptr< CWriteTextFile > CWriteTextFile::MakeWriteTextFile( 
std::unique_ptr< CWriteTextFile > CWriteTextFile::MakeWriteTextFile( 
        const std::string & name, TCompression c )
{
    if( c == COMPRESSION_AUTO ) c = Name2CompressionType( name );

    switch( c ) {
        case COMPRESSION_NONE:
            // return std::auto_ptr< CWriteTextFile >( 
            return std::unique_ptr< CWriteTextFile >( 
                    new CWriteTextFile_CPPStream( name ) );
        default: SRPRISM_ASSERT( false );
    }

    return std::unique_ptr< CWriteTextFile >( nullptr );
    // return std::auto_ptr< CWriteTextFile >( 0 );
}

#undef CHECK_STREAM

END_NS( common )
END_STD_SCOPES

