/*  $Id: tmpstore.cpp 205411 2010-09-17 17:42:11Z morgulis $
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
 * File Description: utilities for temporary file name management
 *
 */

#include <ncbi_pch.hpp>

#include "def.h"

#ifndef WIN32
#	include <sys/types.h>
#	include <unistd.h>

#	define CLOSE close
#	define UNLINK unlink
#	define GETPID getpid
#else
#	include <io.h>
#	include <process.h>

#	define mkstemp(a) _mktemp_s( (a), 7 )
#	define CLOSE _close
#	define UNLINK _unlink
#	define GETPID _getpid
#endif

#include <sstream>
#include <thread>

#include "trace.hpp"
#include "tmpstore.hpp"

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
CTmpStore::CTmpStore( const std::string & tmp_dir_name )
{
    char templ[] = "XXXXXX";
    int fid( -1 );

    if( (fid = mkstemp( templ )) == -1 ) {
        M_TRACE( CTracer::WARNING_LVL, 
                 "could not create unique temporary file suffix" );
    }
    else { CLOSE( fid ); UNLINK( templ ); }

    tmp_name_prefix_ = tmp_dir_name + FPATH_SEP + ".";
    int pid = GETPID();
    auto tid( std::this_thread::get_id() );
    std::ostringstream os;
    os << "." << pid << "." << tid << "." << templ;
    tmp_name_suffix_ = os.str();
}

//------------------------------------------------------------------------------
CTmpStore::~CTmpStore(void)
{
    for( TData::const_iterator i = data_.begin(); i != data_.end(); ++i ) {
        if( UNLINK( CreateName( *i ).c_str() ) < 0 ) {
            M_TRACE( CTracer::WARNING_LVL, "can not unlink " << *i );
        }
    }
}

//------------------------------------------------------------------------------
const std::string CTmpStore::CreateName( const std::string & name )
{ return tmp_name_prefix_ + name + tmp_name_suffix_; }

//------------------------------------------------------------------------------
const std::string CTmpStore::Register( const std::string & name )
{ data_.insert( name ); return CreateName( name ); }

//------------------------------------------------------------------------------
bool CTmpStore::Find( const std::string & name ) const
{ return (data_.find( name ) != data_.end()); }

END_NS( common )
END_STD_SCOPES

