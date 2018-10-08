/*  $Id: memmgr.hpp 351764 2012-02-01 14:07:34Z morgulis $
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

#ifndef __SRPRISM_MEMMGR_HPP__
#define __SRPRISM_MEMMGR_HPP__

#include "../common/def.h"

#include <cstdlib>
#include <vector>
#include <map>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <common/bits.hpp>
#include <srprism/srprismdef.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CMemoryManager
{
    private:

        typedef common::Uint8 TUnit;
        
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode LIMIT     = 1;
            static const TErrorCode ALLOC     = 2;
            static const TErrorCode NOT_FOUND = 3;
            static const TErrorCode REALLOC   = 4;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == LIMIT ) return "memory limit exceeded";
                else if( code == ALLOC ) return "allocation error";
                else if( code == NOT_FOUND ) return "allocated block not found";
                else if( code == REALLOC ) return "reallocation error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException );
        };

        CMemoryManager( TSize memory_limit );
        ~CMemoryManager();

        void * Allocate( TSize request_bytes );
        void Free( void * ptr );
        void * Shrink( void * ptr, TSize request_bytes );

        TSize GetFreeSpace( void ) const { return free_space_*sizeof( TUnit ); }

    private:

        typedef common::Uint8 * TBytePtr;

        // Correspondence of allocated chunks to their sizes (in 8-byte units)
        //
        typedef std::map< TBytePtr, TSize > TAllocMap;

        static const int SHIFT = common::SBinLog< 
            sizeof( TUnit )/sizeof( common::Uint1 ) >::VALUE;

        static TSize Bytes2Units( TSize request )
        { return (request > 0) ? 1 + ((request - 1)>>SHIFT) : 0; }

        TAllocMap::iterator FindBlock( void * ptr );

        TAllocMap alloc_map_;
        TSize free_space_;  // free space in 8-byte words
};

END_NS( srprism )
END_STD_SCOPES

#endif

