/*  $Id: scratch.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: scratch space for bad n-mer filter
 *
 */

#ifndef __SPRISM_SCRATCH_HPP__
#define __SPRISM_SCRATCH_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <cassert>
#include <vector>
#include <algorithm>

#include <common/bits.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <cassert>
#include <vector>
#include <algorithm>

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CScratchBitMap
{
    private:

        typedef common::Uint8 TUnit;
        typedef TUnit * TBuf;

        static const int UNIT_BITS = common::BYTEBITS*sizeof( TUnit );

        static size_t Bits2Units( size_t sz ) { return 1 + (sz - 1)/UNIT_BITS; }

    public:

        CScratchBitMap( void * buf, size_t sz ) 
            : buf_( (TBuf)buf ), buf_sz_( Bits2Units( sz ) )
        { std::fill( buf_, buf_ + buf_sz_, 0 ); }

        void Clear( size_t sz )
        {
            sz = std::min( buf_sz_, Bits2Units( sz ) );
            std::fill( buf_, buf_ + sz, 0 );
        }

        void SetBit( size_t sz, bool set = true )
        {
            SRPRISM_ASSERT( sz < buf_sz_*UNIT_BITS );
            size_t unit( sz/UNIT_BITS ), bit( sz%UNIT_BITS );
            common::AssignBit( bit, buf_[unit], set );
        }

        bool BitIsSet( size_t sz ) const
        { 
            SRPRISM_ASSERT( sz < buf_sz_*UNIT_BITS );
            size_t unit( sz/UNIT_BITS ), bit( sz%UNIT_BITS );
            return common::GetBit( bit, buf_[unit] );
        }

        size_t Size( void ) const { return UNIT_BITS*buf_sz_; }

    private:

         TBuf buf_;
         size_t buf_sz_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

