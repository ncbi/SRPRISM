/*  $Id: seqdef.hpp 537014 2017-05-25 12:32:51Z morgulis $
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
 * File Description: sequence handling utilities
 *
 */

#ifndef __AM_SEQ_SEQDEF_HPP__
#define __AM_SEQ_SEQDEF_HPP__

#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include "../common/bits.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>

#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
/**\brief Types used to symbolically identify sequences. */
typedef std::string TSeqId;

/**\brief Text based sequence descriptions. */
typedef std::string TSeqTitle;

//------------------------------------------------------------------------------
/**\brief Integer type used to enumerate supported encodings. */
typedef common::Uint1 TCoding;

/**\group encodings Symbolic names for encodings. */
/**\@{*/
const TCoding CODING_NCBI2NA = 0;
const TCoding CODING_NCBI4NA = 1;
const TCoding CODING_IUPACNA = 2;
/**\@}*/

//------------------------------------------------------------------------------
/**\brief Integer type used to enumerate strands. */
typedef common::Uint1 TStrand;

/**\brief strands Symbolic strand names. */
/**\@{*/
const TStrand STRAND_FW = 0;
const TStrand STRAND_RV = 1;
/**\@}*/

/**\brief Total number of strands. */
const size_t N_STRANDS = 2;

/**\brief Get a relative (addition modulo 2) strand given 2 strand values.
   
   \param l first strand
   \param r second strand

   \return STRAND_FW is first and second are collinear; STRAND_RV if they
           are opposite to each other.
 */
inline TStrand CombineStrands( TStrand l, TStrand r )
{ return (l == r) ? STRAND_FW : STRAND_RV; }

/**\brief Strand reversal in the usual way. */
inline TStrand ReverseStrand( TStrand s ) { return (TStrand)(1 - s); }

//------------------------------------------------------------------------------
/**\brief Nucleotide base holder type (any encoding). */
typedef common::Uint1 TLetter;

/**\brief Type used to hold sequence sizes and offsets. */
typedef common::Uint4 TSeqSize;

/**\brief Type for nucleotide alphabets represented as STL sequences. */
typedef std::vector< TLetter > TAlphabet;

/**\brief Type for nucleotide alphabets represented as strings. */
typedef std::string TAlphabetString;

/**\brief Type for nucleotide alphabets represented as C arrays. */
typedef const TLetter TAlphabetArray[];

//------------------------------------------------------------------------------
/**\brief Base class for coding properties. */
template< TCoding coding >
struct SCodingTraits_Base {};

/**\group coding_traits Specialization of coding traits. */
/**\@{*/
template<> struct SCodingTraits_Base< CODING_NCBI2NA >
{
    static const size_t LETTER_BITS = 2;
    static TAlphabetArray ALPHABET_ARRAY;
    static TAlphabetArray RC_ALPHABET_ARRAY;
    static TAlphabetArray AMBIG_ARRAY;
};

template<> struct SCodingTraits_Base< CODING_NCBI4NA >
{
    static const size_t LETTER_BITS = 4;
    static TAlphabetArray ALPHABET_ARRAY;
    static TAlphabetArray RC_ALPHABET_ARRAY;
    static TAlphabetArray AMBIG_ARRAY;
    static const TLetter GAP_LETTER = '\000';
};

template<> struct SCodingTraits_Base< CODING_IUPACNA >
{
    static const size_t LETTER_BITS = 8;
    static TAlphabetString ALPHABET_STRING;
    static TAlphabetString RC_ALPHABET_STRING;
    static TAlphabetString AMBIG_STRING;
    static const TLetter GAP_LETTER = '-';
};
/**\@}*/

/**\brief Derived coding properties and properties that do not have
          to be specialized. 
 */
template< TCoding coding >
struct SCodingTraits : public SCodingTraits_Base< coding >
{
    private: typedef SCodingTraits_Base< coding > TBase;

    public:

        static const unsigned int PACK_FACTOR = 
            common::BYTEBITS/TBase::LETTER_BITS;

        static const size_t LETTER_SHIFT = 
            common::SBinLog< TBase::LETTER_BITS >::VALUE;

        static const size_t NUM_VALUES = 
            common::SBitFieldTraits< size_t, TBase::LETTER_BITS >::MAX + 1;

        static TAlphabet       ALPHABET;
        static bool            NAMBIG[SCodingTraits< coding >::NUM_VALUES];
        static TLetter         RC[SCodingTraits< coding >::NUM_VALUES];

        /**\brief Convert letters between encodings. */
		template< TCoding dst_code > static TLetter Recode( TLetter l );
};

/**\brief Byte-level reverse complement table. */
extern common::Uint1 
BYTE_RC[common::SBitFieldTraits< size_t, common::BYTEBITS >::MAX + 1];

/**\brief Initialization of various coding tables.

   This function needs to be called before any other facilities 
   for work with encodings are used.
*/
void InitCoding(void);

//------------------------------------------------------------------------------
struct SSeqData_Base { 
    SSeqData_Base( TSeqSize size ) : size( size ) {}

    TSeqSize size;
};

template< typename word_t >
struct SSeqData_TBase : public SSeqData_Base
{
    typedef word_t TWord;
    typedef std::vector< TWord > TSeq;

    SSeqData_TBase( TSeq & seq, TSeqSize size ) 
        : SSeqData_Base( size ), seq( seq ) {}

    TSeq & seq;
};

template< TCoding coding, typename word_t >
struct SSeqData : public SSeqData_TBase< word_t >
{
    typedef SSeqData_TBase< word_t > TBase;
    typedef typename TBase::TWord TWord;
    typedef typename TBase::TSeq TSeq;

    static const TCoding CODING = coding;
    typedef SCodingTraits< CODING > TCodingTraits;

    SSeqData( TSeq & seq, TSeqSize size ) : TBase( seq, size ) {}

    size_t Ambig( common::Sint2 start, common::Sint2 end ) const
    {
        size_t result( 0 );
        TSeq & seq( this->seq );

        for( TSeqSize i( std::max( 0, (int)start) ),
                      sz( std::min( this->size, (TSeqSize)end ) ); 
             i < sz; ++i ) 
        {
            if( !TCodingTraits::NAMBIG[seq[i]] ) ++result;
        }

        return result;
    }
};

//------------------------------------------------------------------------------
template< TSeqSize S, TSeqSize E, TCoding coding, typename unit_t >
inline unit_t SelectLetters( unit_t u )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t EBIT  = UBITS - S*LBITS;

    common::SelectBits< EBIT - (E-S)*LBITS, EBIT >( u );
    return u;
}

template< TCoding coding, typename unit_t >
inline unit_t SelectLetters( TSeqSize s, TSeqSize e, unit_t u )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    if( s >= e ) return 0;
    size_t ebit( UBITS - (s<<LSHIFT) );
    common::SelectBits( u, ebit - ((e - s)<<LSHIFT), ebit );
    return u;
}

//------------------------------------------------------------------------------
template< TSeqSize S, TSeqSize E, TCoding coding, typename unit_t >
inline void AssignLetters_Unsafe( unit_t & du, unit_t su )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t EBIT  = UBITS - S*LBITS;
    static const size_t SBIT  = EBIT - (E-S)*LBITS;
    SRPRISM_ASSERT( (su&common::RMask< SBIT, EBIT >( su ) == 0) );
    common::AssignBits_Unsafe< SBIT, EBIT >( du, su );
}

template< TCoding coding, typename unit_t >
inline void AssignLetters_Unsafe( 
        TSeqSize s, TSeqSize e, unit_t & du, unit_t su )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    size_t ebit( UBITS - (s<<LSHIFT) );
    size_t sbit( ebit - ((e - s)<<LSHIFT) );
    SRPRISM_ASSERT( ((su&(common::RMask< unit_t >( sbit, ebit ))) == 0) );
    common::AssignBits_Unsafe( sbit, ebit, du, su );
}

template< TSeqSize S, TSeqSize E, TCoding coding, typename unit_t >
inline void AssignLetters( unit_t & du, unit_t su )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t EBIT  = UBITS - S*LBITS;
    static const size_t SBIT  = EBIT - (E-S)*LBITS;
    common::AssignBits< SBIT, EBIT >( du, su );
}

template< TCoding coding, typename unit_t >
inline void AssignLetters( TSeqSize s, TSeqSize e, unit_t & du, unit_t su )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    size_t ebit( UBITS - (s<<LSHIFT) );
    size_t sbit( ebit - ((e - s)<<LSHIFT) );
    common::AssignBits( sbit, ebit, du, su );
}

//------------------------------------------------------------------------------
template< TCoding coding, typename unit_t >
inline unit_t GetLetters_LS( TSeqSize s, TSeqSize e, unit_t u )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    if( s >= e ) return 0;
    size_t shift( UBITS - (e<<LSHIFT) );
    u = (SelectLetters< coding >( s, e, u )>>shift);
    return u;
}

//------------------------------------------------------------------------------
template< TCoding coding, typename unit_t >
inline TLetter GetLetter( unit_t u, TSeqSize i )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    size_t ebit( UBITS - (i<<LSHIFT) );
    return (TLetter)common::GetField( ebit - LBITS, ebit, u );
}

template< TCoding coding, typename iter_t >
inline TLetter GetStreamLetter( iter_t i, TSeqSize pos )
{
    typedef typename common::SIterTraits< iter_t >::TValue TValue;
    typedef SCodingTraits< coding > TTraits;
    static const TSeqSize ULETTERS = sizeof( TValue )*TTraits::PACK_FACTOR;
    static const size_t SHIFT = common::SBinLog< ULETTERS >::VALUE;
    static const TValue MASK = common::SBitFieldTraits< TValue, SHIFT >::MASK;

    TValue v( *(i + (pos>>SHIFT)) ); pos &= MASK;
    return GetLetter< coding >( v, pos );
}

//------------------------------------------------------------------------------
template < TCoding coding, typename unit_t >
inline void SetLetter( unit_t & u, TSeqSize i, TLetter l )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t UBITS = common::SIntTraits< unit_t >::BITS;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

    SRPRISM_ASSERT( (l>>LBITS) == 0 );
    size_t ebit( UBITS - (i<<LSHIFT) );
    common::SetField_Unsafe( ebit - LBITS, ebit, u, (unit_t)l );
}

template< TCoding coding, typename iter_t >
inline void SetStreamLetter( iter_t i, TSeqSize pos, TLetter l )
{
    typedef typename common::SIterTraits< iter_t >::TValue TValue;
    typedef SCodingTraits< coding > TTraits;
    static const TSeqSize ULETTERS = sizeof( TValue )*TTraits::PACK_FACTOR;
    static const size_t SHIFT = common::SBinLog< ULETTERS >::VALUE;
    static const TValue MASK = common::SBitFieldTraits< TValue, SHIFT >::MASK;

    TValue & v( *(i + (pos>>SHIFT)) ); pos &= MASK;
    SetLetter< coding >( v, pos, l );
}

//------------------------------------------------------------------------------
template< TCoding coding, size_t NLETTERS, typename unit_t >
inline void PushLetter( unit_t & u, TLetter l )
{
    typedef SCodingTraits< coding > Traits;
    static const size_t LBITS = Traits::LETTER_BITS;
    static const size_t NBITS = LBITS*NLETTERS;
    static const unit_t MASK = common::SBitFieldTraits< unit_t, NBITS >::MASK;

    common::PushBits_Unsafe< LBITS >( u, (unit_t)l );
    u &= MASK;
}

//------------------------------------------------------------------------------
template< 
    TCoding d_coding, TCoding s_coding, 
    typename d_iter_t, typename s_iter_t >
inline void Recode( d_iter_t di, s_iter_t si, TSeqSize len )
{
    for( TSeqSize i = 0; i < len; ++i ) {
        TLetter l( GetStreamLetter< s_coding >( si, i ) );
        SetStreamLetter< d_coding >( 
                di, i, 
                SCodingTraits< s_coding >::template Recode< d_coding >( l ) );
    }
}

//------------------------------------------------------------------------------
template< TCoding coding >
struct SReverseComplement
{
    typedef SCodingTraits< coding > TTraits;

    template< typename unit_t >
    void operator()( unit_t & d, unit_t s, TSeqSize len )
    {
        static const TSeqSize ULETTERS = sizeof( unit_t )*TTraits::PACK_FACTOR;

        SRPRISM_ASSERT( len < ULETTERS );

        for( TSeqSize i = 0; i < len; ++i ) {
            TLetter l( GetLetter< coding >( s, i ) );
            SetLetter< coding >( d, len - i - 1, TTraits::RC[l] );
        }
    }
};

template<>
struct SReverseComplement< CODING_NCBI2NA >
{
    static const TCoding CODING = CODING_NCBI2NA;
    typedef SCodingTraits< CODING > TTraits;

    template< typename unit_t >
    void operator()( unit_t & d, unit_t s )
    {
        static const size_t UBITS = common::BYTEBITS*sizeof( unit_t );
        static const size_t MASK = 
            common::SBitFieldTraits< unit_t, common::BYTEBITS >::MASK;

        for( size_t i = 0; i < UBITS; i += common::BYTEBITS ) {
            common::Uint1 sbyte( (common::Uint1)((s>>i)&MASK) );
            common::PushBits_Unsafe< common::BYTEBITS >( 
                    d, (unit_t)BYTE_RC[sbyte] );
        }
    }

    template< typename unit_t >
    void operator()( unit_t & d, unit_t s, TSeqSize len )
    {
        static const size_t UBYTES = sizeof( unit_t );
        static const size_t UBITS = common::BYTEBITS*UBYTES;

#ifndef NDEBUG
        static const TSeqSize ULETTERS = sizeof( unit_t )*TTraits::PACK_FACTOR;
#endif

        static const size_t LSHIFT = 
            common::SBinLog< TTraits::LETTER_BITS >::VALUE;

#ifndef NDEBUG
        SRPRISM_ASSERT( len <= ULETTERS );
#endif

        size_t shift( UBITS - (len<<LSHIFT) );
        (*this)( d, s );
        d <<= shift;
    }

    template< typename unit_t >
    void operator()( 
            unit_t * d, const unit_t * s, TSeqSize start, TSeqSize len )
    {
        static const size_t UBYTES = sizeof( unit_t );
        static const size_t UBITS = common::BYTEBITS*UBYTES;
        static const size_t LSHIFT = 
            common::SBinLog< TTraits::LETTER_BITS >::VALUE;
        static const TSeqSize ULETTERS = sizeof( unit_t )*TTraits::PACK_FACTOR;
        static const size_t LLSHIFT = common::SBinLog< ULETTERS >::VALUE;
        static const TSeqSize LLMASK = 
            common::SBitFieldTraits< TSeqSize, LLSHIFT >::MASK;

        TSeqSize n_units( (len + ULETTERS - 1)/ULETTERS );
        s += (start>>LLSHIFT);
        len += (start&LLMASK);
        TSeqSize lenrem( len&LLMASK );
        const unit_t * se( s + (len>>LLSHIFT) );

        if( lenrem > 0 ) {
            TSeqSize lshift( lenrem<<LSHIFT ),
                     rshift( UBITS - lshift );

            while( se > s ) {
                unit_t su( ((*se)>>rshift) + ((*(se-1))<<lshift) );
                operator()( *d++, su );
                --se;
                --n_units;
            }

            if( n_units > 0 ) operator()( *d, *se, lenrem );
        }
        else while( se > s ) operator()( *d++, *--se );
    }
};

template< TCoding coding, typename unit_t >
void ReverseComplement( unit_t & d, unit_t s, TSeqSize len )
{ SReverseComplement< coding >()( d, s, len ); }

template< TCoding coding, typename unit_t >
void ReverseComplement( unit_t & d, unit_t s )
{ SReverseComplement< coding >()( d, s ); }

template< TCoding coding, typename unit_t >
void ReverseComplement( 
        unit_t * d, const unit_t * s, TSeqSize start, TSeqSize len )
{ SReverseComplement< coding >()( d, s, start, len ); }

//------------------------------------------------------------------------------
template< TCoding coding, typename word_t >
word_t GetWord( const word_t * src, TSeqSize off )
{
    typedef SCodingTraits< coding > TTraits;
    static const size_t BITSHIFT = 
        common::SBinLog< TTraits::LETTER_BITS >::VALUE;
    static const size_t WORD_BITS = sizeof( word_t )*common::BYTEBITS;
    static const size_t WSHIFT = common::SBinLog< WORD_BITS >::VALUE;
    static const word_t WMASK = 
        common::SBitFieldTraits< word_t, WSHIFT >::MASK;

    off <<= BITSHIFT;
    src += (off>>WSHIFT);
    off &= WMASK;
    word_t res( (*src)<<off );
    return (off == 0) ? res : res + ((*(src+1))>>(WORD_BITS - off));
}

//------------------------------------------------------------------------------
template< typename word_t >
size_t NTriplets( word_t word, size_t len, common::Uint8 & triplets_present )
{
    size_t res( 0 );

    if( len > 2 ) {
        for( int i = 0; (size_t)i < len - 2; ++i ) {
            word_t triplet( 
                    GetLetters_LS< CODING_NCBI2NA >( i, i + 3, word ) );

            if( !common::GetBit( triplet, triplets_present ) ) {
                ++res;
                common::SetBit( triplets_present, triplet );
            }
        }
    }

    return res;
}

extern const std::string CS_ALPHABET_STRING;
extern std::string Color2IUPACNA( const std::string & c_str, TLetter lead );
extern std::string IUPACNA2Color( const std::string & str, TLetter lead );

END_NS( seq )
END_STD_SCOPES

#endif

