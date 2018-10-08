/*  $Id: seqiter.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: sequence iterator over a packed sequence buffer
 *
 */

#ifndef __SRPRISM_SEQITER_HPP__
#define __SRPRISM_SEQITER_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <utility>
#include <algorithm>

#include <common/bits.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <utility>
#include <algorithm>

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
template< typename word_t, typename unit_t >
class CBitIterator
{
    public:

        typedef word_t TWord;
        typedef unit_t TUnit;

    private:

        static const size_t WBITS = common::BYTEBITS*sizeof( TWord );
        static const size_t UBITS = common::BYTEBITS*sizeof( TUnit );
        static const size_t WSHIFT = common::SBinLog< WBITS >::VALUE;
        static const size_t WMASK = 
            common::SBitFieldTraits< size_t, WSHIFT >::MASK;

    public:

        CBitIterator( const TWord * data, size_t off, size_t len )
            : data_( data ), n_words_( 1 + ((len - 1)>>WSHIFT) ),
              off_( off ), len_( len )
        {
            SRPRISM_ASSERT( data_ != 0 );
            SRPRISM_ASSERT( off <= len );
            SRPRISM_ASSERT( sizeof( TWord ) >= sizeof( TUnit ) );
        }

    private:

        void Advance( size_t adj )
        { off_ = std::min( off_ + adj, len_ ); }

        void Retreat( size_t adj )
        { off_ -= std::min( adj, off_ ); }

    public:

        bool End( void ) const { return (off_ == len_); }
        bool Begin( void ) const { return (off_ == 0 ); }

        CBitIterator & operator+=( ssize_t adj )
        {
            if( adj > 0 ) Advance( adj ); else Retreat( -adj );
            return *this;
        }

        CBitIterator & operator-=( ssize_t adj )
        { return operator+=( -adj ); }

        CBitIterator & operator++( void ) { return operator+=( 1 ); }
        CBitIterator & operator--( void ) { return operator-=( 1 ); }

        std::pair< TUnit, size_t > operator*( void ) const
        {
            size_t l( std::min( len_ - off_, (size_t)UBITS ) );
            size_t w( off_>>WSHIFT ), b( off_&WMASK );
            const TWord & word( data_[w] );
            
            if( b + UBITS > WBITS ) {
                TUnit res( (TUnit)common::GetFieldBack( 
                            WBITS - b, word )<<(UBITS - (WBITS - b)) );

                if( w + 1 >= n_words_ ) return std::make_pair( res, l );
                else {
                    const TWord & word1( data_[w+1] );
                    TUnit res1( (TUnit)common::GetField( 
                                WBITS + WBITS - UBITS - b, WBITS, word1 ) );
                    return std::make_pair( res + res1, l );
                }
            }
            else {
                TUnit res( (TUnit)(word>>(WBITS - UBITS - b)) );
                return std::make_pair( res, l );
            }
        }

    protected:

        const TWord * data_;
        size_t n_words_;
        size_t off_;
        size_t len_;
};

//------------------------------------------------------------------------------
template< TCoding CODING, typename word_t, typename unit_t > 
class CSeqRvIterator;

template< TCoding CODING, typename word_t, typename unit_t >
class CSeqFwIterator : public CBitIterator< word_t, unit_t >
{
    private:

        typedef CBitIterator< word_t, unit_t > TBase;
        static const size_t LBITS = seq::SCodingTraits< CODING >::LETTER_BITS;
        static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

        typedef typename TBase::TWord TWord;

    public:

        typedef CSeqRvIterator< CODING, word_t, unit_t > TReverse;

        CSeqFwIterator( const TWord * data, TSeqSize off, TSeqSize len )
            : TBase( data, (size_t)(off<<LSHIFT), (size_t)(len<<LSHIFT) )
        {}

        CSeqFwIterator( const TBase & base_it ) : TBase( base_it ) {}

        CSeqFwIterator & operator+=( ssize_t adj )
        { TBase::operator+=( adj<<LSHIFT ); return *this; }

        CSeqFwIterator & operator-=( ssize_t adj )
        { return operator+=( -adj ); }

        CSeqFwIterator & operator++( void ) { return operator+=( 1 ); }
        CSeqFwIterator & operator--( void ) { return operator-=( 1 ); }

        std::pair< unit_t, size_t > operator*( void ) const
        {
            std::pair< unit_t, size_t > r( TBase::operator*() );
            r.second >>= LSHIFT;
            return r;
        }
};

//------------------------------------------------------------------------------
template< TCoding CODING, typename word_t, typename unit_t >
class CSeqRvIterator : public CBitIterator< word_t, unit_t >
{
    private:

        typedef CBitIterator< word_t, unit_t > TBase;

    public:

        typedef typename TBase::TWord TWord;
        typedef typename TBase::TUnit TUnit;

    private:

        static const size_t LBITS = seq::SCodingTraits< CODING >::LETTER_BITS;
        static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;
        static const size_t ULETTERS = 
            sizeof( TUnit )*seq::SCodingTraits< CODING >::PACK_FACTOR;
        static const size_t UBITS = common::BYTEBITS*sizeof( TUnit );
        static const size_t WBITS = common::BYTEBITS*sizeof( TWord );

    public:

        typedef CSeqFwIterator< CODING, word_t, unit_t > TReverse;

        CSeqRvIterator( const TWord * data, TSeqSize off, TSeqSize len )
            : TBase( data, (size_t)(off<<LSHIFT), (size_t)(len<<LSHIFT) )
        {}

        CSeqRvIterator( const TBase & base_it ) : TBase( base_it ) {}

        bool End( void ) const { return TBase::Begin(); }
        bool Begin( void ) const { return TBase::End(); }

        CSeqRvIterator & operator+=( ssize_t adj )
        { TBase::operator-=( adj<<LSHIFT ); return *this; }

        CSeqRvIterator & operator-=( ssize_t adj )
        { return operator+=( -adj ); }

        CSeqRvIterator & operator++( void ) { return operator+=( 1 ); }
        CSeqRvIterator & operator--( void ) { return operator-=( 1 ); }

        std::pair< TUnit, size_t > operator*( void ) const
        {
            if( this->off_ >= UBITS ) {
                TBase t( *this ); t -= UBITS;
                TUnit res( (*t).first ), res_rc( 0 );
                seq::ReverseComplement< CODING >( res_rc, res );
                return std::make_pair( res_rc, ULETTERS );
            }
            else {
                TUnit res( (TUnit)common::GetField( 
                            WBITS - UBITS, WBITS, *(this->data_) ) ),
                      res_rc( 0 );
                seq::ReverseComplement< CODING >( res_rc, res );
                return std::make_pair( 
                        (res_rc<<(UBITS - this->off_)), (this->off_>>LSHIFT) );
            }
        }
};

END_NS( srprism )
END_STD_SCOPES

#endif

