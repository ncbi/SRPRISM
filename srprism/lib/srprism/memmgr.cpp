/*  $Id: memmgr.cpp 426067 2014-02-05 18:25:10Z morgulis $
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
 * File Description: memory management code
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <cassert>
#include <cstdlib>
#include <algorithm>

#include "../common/trace.hpp"
#include "memmgr.hpp"

//
// uncomment next line to enable allocation tracing
// #define M_TRACE_ALLOC_ENABLED 1

#ifdef M_TRACE_ALLOC_ENABLED
#   define M_TRACE_ALLOC M_TRACE
#else
#   define M_TRACE_ALLOC(lvl,msg)
#endif

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
CMemoryManager::CMemoryManager( TSize memory_limit )
    : free_space_( Bytes2Units( memory_limit ) )
{
    SRPRISM_ASSERT( free_space_ > 0 );
    M_TRACE_ALLOC( CTracer::INFO_LVL, "CMemoryManager(): " << free_space_ );
}

//------------------------------------------------------------------------------
CMemoryManager::~CMemoryManager()
{
    while( !alloc_map_.empty() ) {
        TBytePtr p( alloc_map_.begin()->first );
        alloc_map_.erase( alloc_map_.begin() );
        free( p );
    }

    M_TRACE_ALLOC( CTracer::INFO_LVL, "~CMemoryManager()" );
}

//------------------------------------------------------------------------------
void * CMemoryManager::Allocate( TSize request_bytes )
{
    if( request_bytes == 0 ) return 0;
    TSize units( Bytes2Units( request_bytes ) );

    if( units > free_space_ ) {
        M_THROW( CException, LIMIT, 
                 "request: " << request_bytes << " bytes; limit: " <<
                 free_space_*sizeof( TUnit ) << " bytes" );
    }

    request_bytes = units*sizeof( TUnit );
    void * result( malloc( request_bytes ) );

    if( result == 0 ) {
        M_THROW( CException, ALLOC, "request: " << request_bytes << " bytes" );
    }

    alloc_map_[(TBytePtr)result] = units;
    free_space_ -= units;
    M_TRACE_ALLOC( CTracer::INFO_LVL, 
             "Allocate(): " << units << " (" << free_space_ << " free)" );
    return result;
}

//------------------------------------------------------------------------------
inline CMemoryManager::TAllocMap::iterator 
CMemoryManager::FindBlock( void * ptr )
{
    TAllocMap::iterator i( alloc_map_.find( (TBytePtr)ptr ) );

    if( i == alloc_map_.end() ) {
        M_THROW( CException, NOT_FOUND, 
                 "address: " << std::hex << (Uint8)ptr );
    }

    return i;
}

//------------------------------------------------------------------------------
void CMemoryManager::Free( void * ptr )
{
    if( ptr == 0 ) return;
    TAllocMap::iterator i( FindBlock( ptr ) );
    TSize units( i->second );
    alloc_map_.erase( i );
    free( ptr );
    free_space_ += units;
    M_TRACE_ALLOC( CTracer::INFO_LVL, 
             "De-Allocate(): " << units << " (" << free_space_ << " free)" );
}

//------------------------------------------------------------------------------
void * CMemoryManager::Shrink( void * ptr, TSize request_bytes )
{
    if( ptr == 0 ) return 0;
    if( request_bytes == 0 ) { Free( ptr ); return 0; }
    TSize new_units( Bytes2Units( request_bytes ) );
    TAllocMap::iterator i( FindBlock( ptr ) );
    TSize old_units( i->second );
    SRPRISM_ASSERT( new_units <= old_units );
    if( new_units == old_units ) return ptr;
    request_bytes = new_units*sizeof( TUnit );
    ptr = realloc( ptr, request_bytes );

    if( ptr == 0 ) {
        M_THROW( CException, REALLOC, 
                 "address: " << std::hex << (Uint8)i->first << std::dec <<
                 "; request: " << request_bytes << " bytes" );
    }

    alloc_map_.erase( i );
    alloc_map_[(TBytePtr)ptr] = new_units;
    free_space_ += (old_units - new_units);
    M_TRACE_ALLOC( CTracer::INFO_LVL, 
             "Re-Allocate(): " << old_units << " ---> " << new_units <<
             " (" << free_space_ << " free)" );
    return ptr;
}

END_NS( srprism )
END_STD_SCOPES
