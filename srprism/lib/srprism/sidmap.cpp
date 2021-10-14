/*  $Id: sidmap.cpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: subject sequence id map interface
 *
 */

#include <ncbi_pch.hpp>

#include <cstdlib>
#include <algorithm>

#include "../common/trace.hpp"
#include "../common/textfile.hpp"
#include "sidmap.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
const char * CSIdMap::FILE_SFX = ".imp";

//------------------------------------------------------------------------------
void CSIdMap::CleanUp(void)
{
    if( data_ != 0 ) { mem_mgr_.Free( data_ ); data_ = 0; }
    if( offs_ != 0 ) { mem_mgr_.Free( offs_ ); offs_ = 0; }
}

//------------------------------------------------------------------------------
CSIdMap::CSIdMap( const std::string & name, CMemoryManager & mem_mgr )
    : mem_mgr_( mem_mgr ), offs_( 0 ), data_( 0 ), n_ids_( 0 )
{
    try{
        std::string fname( name + FILE_SFX );
        // std::auto_ptr< CReadTextFile > fidmap( 
        std::unique_ptr< CReadTextFile > fidmap( 
                common::CReadTextFile::MakeReadTextFile( fname ) );

        if( fidmap->Eof() ) {
            M_THROW( CException, FORMAT, "empty id map " << fname );
        }

        n_ids_ = atol( fidmap->GetLine().c_str() );
        M_TRACE( CTracer::INFO_LVL,
                 "loading ids for " << n_ids_ << " sequences" );

        offs_ = (size_t *)mem_mgr_.Allocate( sizeof( size_t )*(n_ids_ + 1) );
        std::fill( offs_, offs_ + n_ids_ + 1, 0 );
        size_t data_size = mem_mgr_.GetFreeSpaceSize();
        data_ = (char *)mem_mgr_.Allocate( mem_mgr_.GetFreeSpaceSize() );
        size_t total( 0 );

        for( size_t lc( 0 ); lc < n_ids_; ++lc ) {
            if( fidmap->Eof() ) {
                M_THROW( CException, FORMAT,
                         "end of file reached before reading requested "
                         "number " << n_ids_ << " of ids; at line " <<
                         fidmap->LineNo() );
            }

            std::string id( fidmap->GetLine() );

            if( total + id.size() > data_size ) {
                M_THROW( CException, MEMORY, "at line " << fidmap->LineNo() );
            }

            offs_[lc] = total;
            std::copy( id.begin(), id.end(), data_ + total );
            total += id.size();
        }

        offs_[n_ids_] = total;
        data_ = (char *)mem_mgr_.Shrink( data_, total );
    }
    catch( std::exception & e ) {
        CleanUp();
        M_TRACE( CTracer::ERROR_LVL,
                 "error loading database id map; "
                 "subject ids will not appear in the output [" <<
                 e.what() << "]" );
    }
}

//------------------------------------------------------------------------------
CSIdMap::~CSIdMap(void)
{ CleanUp(); }

END_NS( srprism )
END_STD_SCOPES

