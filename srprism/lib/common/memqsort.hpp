/*  $Id: memqsort.hpp 187432 2010-03-31 16:11:09Z morgulis $
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
 * File Description: quicksort implementation
 *
 */

#ifndef __AM_COMMON_MEMQSORT_HPP__
#define __AM_COMMON_MEMQSORT_HPP__

#include <common/def.h>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
template< typename ptr_t, typename comp_t, typename swap_t >
inline ptr_t MemQPart( 
        ptr_t start, ptr_t end, size_t size, 
        const comp_t & less, const swap_t & swapper )
{
    ptr_t i = start + size, j = end - size;

    while( i < j ) {
        while( i < j && !less( j, start, size ) ) j -= size;
        while( i < j && less( i, start, size ) ) i += size;
        if( i < j ) { swapper( i, j, size ); i += size; j -= size; }
    }

    if( less( i, start, size ) ) swapper( i, start, size );
    else if( start + size != i ) swapper( i - size, start, size );
    return i;
}

//------------------------------------------------------------------------------
template< typename ptr_t, typename comp_t, typename swap_t >
void MemQSort( 
        ptr_t start, ptr_t end, size_t size, 
        const comp_t & less, const swap_t & swapper )
{
    if( start + size < end ) {
        ptr_t i = MemQPart( start, end, size, less, swapper );
        MemQSort( start, i, size, less, swapper );
        MemQSort( i, end, size, less, swapper );
    }
}

END_NS( common )
END_STD_SCOPES

#endif

