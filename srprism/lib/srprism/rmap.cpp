/*  $Id: rmap.cpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: repeat map interface
 *
 */

#include <ncbi_pch.hpp>

#include <algorithm>

#include "../common/binfile.hpp"
#include "../common/trace.hpp"
#include "rmap.hpp"

START_STD_SCOPES
START_NS( srprism )

USE_NS( common )

//------------------------------------------------------------------------------
//
// for backwords compatibility, if we can not read the file we just keep
// the structure empty and give a warning
//
// the repeat map file is a sequence of binary records of length 8; integers
// are stored in machine-dependent byte-order; each record has the following
// structure:
//
// bytes 0-3:   16-mer value;
// bytes 4-5;   binary log of repeat count for forward strand;
// bytes 6-7:   binary log of repeat count for reverse strand;
//
CRMap::CRMap( const std::string & name )
{
    std::string repname( name + CIndexBase::IDX_REPMAP_SFX );
    M_TRACE( CTracer::INFO_LVL, "reading repeat counts from " << repname );

    try {
        CReadBinFile repmap( repname );

        while( !repmap.Eof() ) {
            TPrefix p( 0 ); Uint2 f( 0 ), r( 0 );
            if( repmap.Read( (char *)&p, sizeof( TPrefix ), true ) == 0 ) break;
            repmap.Read( (char *)&f, sizeof( Uint2 ), true );
            repmap.Read( (char *)&r, sizeof( Uint2 ), true );
            SRPRISM_ASSERT( prefixes_.empty() || p > *prefixes_.rbegin() );
            prefixes_.push_back( p );
            ranks_.push_back( (Uint1)std::max( f, r ) );
        }
    }
    catch( CFileBase::CException & e ) {
        prefixes_.clear();
        ranks_.clear();
        M_TRACE( CTracer::WARNING_LVL,
                 "could not import the repeat map file: " << e.what() );
    }

    M_TRACE( CTracer::INFO_LVL, 
             "got repeat counts for " << prefixes_.size() << " 16-mers" );
}

END_NS( srprism )
END_STD_SCOPES

