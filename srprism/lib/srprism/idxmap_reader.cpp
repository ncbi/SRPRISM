/*  $Id: idxmap_reader.cpp 205423 2010-09-17 19:16:42Z morgulis $
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

#include <ncbi_pch.hpp>

#include "../common/trace.hpp"
#include "idxmap_reader.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
CIdxMapReader::CIdxMapReader( const std::string & name )
    : buf_( BUFSIZE, 0 ), is_( name ), sz_( 0 ), start_( 0 ), curr_( 0 )
{
    M_TRACE( CTracer::INFO_LVL, "opened index map " << name );
}

//------------------------------------------------------------------------------
inline bool CIdxMapReader::NextBuf(void)
{
    start_ += sz_;
    CReadBinFile::TSize bytes = 
        is_.Read( (char *)&buf_[0], BUFSIZE*sizeof( TOffset ) );
    sz_ = bytes/sizeof( TOffset ); curr_ = 0;
    
    if( bytes == 0 && start_ != NUM_UNITS ) {
        M_THROW( CException, SIZE, "unexpected end of index map" );
    }

    return (bytes != 0);
}

//------------------------------------------------------------------------------
bool CIdxMapReader::Seek( TUnit target )
{
    while( true ) {
        if( curr_ == sz_ ) { if( !NextBuf() ) return false; }
        else if( start_ + curr_ < target || buf_[curr_] == 0 ) ++curr_;
        else if( buf_[curr_] == END_VAL ) return false;
        else return true;
    }
}

END_NS( srprism )
END_STD_SCOPES

