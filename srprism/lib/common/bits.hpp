/*  $Id: bits.hpp 254996 2011-02-18 16:36:04Z morgulis $
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
 * File Description: utilities for bit operations
 *
 */

#ifndef __AM_COMMON_BITS_HPP__
#define __AM_COMMON_BITS_HPP__

#include "../common/def.h"

/*
#ifndef NCBI_CPP_TK

#include <common/def.h>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#endif
*/

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
/**\brief Compile time computation of binary logarithm. */
template< size_t V > struct SBinLog;

/**\brief Binary log is not defined for argument 0. */
template<> struct SBinLog< 0 > {};

/**\brief Binary log (1) = 0 as induction base. */
template<> struct SBinLog< 1 > { static const size_t VALUE = 0; };

/**\brief Recursive binary log computation for arguments > 1. */
template< size_t V > struct SBinLog
{
    static const size_t VALUE = 1 + SBinLog< (V>>1) >::VALUE;
};

inline size_t BinLog( size_t v )
{
    if( v == 1 ) return 0;
    else return 1 + BinLog( v>>1 );
}

//------------------------------------------------------------------------------
/**\brief Properties of integer values fitting into constant number of bits. */
template< typename int_t, size_t bits >
struct SBitFieldTraits
{
    public:

        /**\brief Base unsigned integer type. */
        typedef int_t TType;

    private:

        /**\brief Maximum supported template parameter value. */
        static const size_t MAX_BITS = BYTEBITS*sizeof( TType );

    public:

        /**\brief Maximum value of the bitfield type. */
        static const TType MAX = (((TType)~0ULL)>>(MAX_BITS - bits));
    
        /**\brief Mask needed to trim a value to the bitfield type. */
        static const TType MASK = MAX;

        /**\brief Bit flipped mask. */
        static const TType RMASK = (TType)~MASK;
};

template< typename int_t >
struct SBitFieldTraits< int_t, 0 >
{
	public: 
		
		typedef int_t TType;
		static const TType MAX = 0;
		static const TType MASK = 0;
		static const TType RMASK = (TType)~MASK;
};

//------------------------------------------------------------------------------
/**\brief Properties of unsigned integer types. */
template< typename int_t >
struct SIntTraits
{
    /**\brief Base unsigned integer type. */
    typedef int_t TType;

    /**\brief Number of bits used by values of the type. */
    static const size_t BITS = BYTEBITS*sizeof( TType );

    /**\brief Maximum value of the type. */
    static const TType MAX = SBitFieldTraits< TType, BITS >::MAX;
};

//------------------------------------------------------------------------------
template< typename iter_t >
struct SIterTraits
{
    typedef iter_t TIter;
    typedef typename TIter::TValue TValue;
};

template< typename type_t >
struct SIterTraits< type_t * >
{
    typedef type_t * TIter;
    typedef type_t TValue;
};

//------------------------------------------------------------------------------
/**\brief Properties of bit ranges within unsigned integer values. */
template< typename int_t, size_t SBIT, size_t EBIT >
struct SBitRangeTraits
{
    /**\brief Base unsigned integer type. */
    typedef int_t TType;

    /**\brief Bit mask corresponding to the bit range. */
    static const TType MASK = 
        ((SBitFieldTraits< TType, EBIT - SBIT >::MASK)<<SBIT);

    /**\brief Reverse mask. */
    static const TType RMASK = (TType)~MASK;
};

//------------------------------------------------------------------------------
/**\brief Properties of a single bit within an unsigned integer value. */
template< typename int_t, size_t BIT >
struct SBitTraits : public SBitRangeTraits< int_t, BIT, BIT + 1 > {};

//------------------------------------------------------------------------------
template< typename int_t >
inline int_t MaskBack( size_t bits )
{ 
    return (bits == SIntTraits< int_t >::BITS) ? ~(int_t)0 
                                               : ((int_t)1<<bits) - 1; 
}

template< typename int_t >
inline int_t RMaskBack( size_t bits )
{ return ~MaskBack< int_t >( bits ); }

template< typename int_t >
inline int_t RMaskFront( size_t bits )
{ return MaskBack< int_t >( SIntTraits< int_t >::BITS - bits ); }

template< typename int_t >
inline int_t MaskFront( size_t bits )
{ return ~RMaskFront< int_t >( bits ); }

template< typename int_t >
inline int_t Mask( size_t sbit, size_t ebit )
{ return (MaskBack< int_t >( ebit - sbit )<<sbit); }

template< typename int_t >
inline int_t RMask( size_t sbit, size_t ebit )
{ return ~Mask< int_t >( sbit, ebit ); }

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline void SetBits( int_t & v )
{ v |= SBitRangeTraits< int_t, SBIT, EBIT >::MASK; }

template< typename int_t >
inline void SetBitsBack( int_t & v, size_t bits )
{ v |= MaskBack< int_t >( bits ); }

template< typename int_t >
inline void SetBitsFront( int_t & v, size_t bits )
{ v |= MaskFront< int_t >( bits ); }

template< typename int_t >
inline void SetBits( int_t & v, size_t sbit, size_t ebit )
{ v |= Mask< int_t >( sbit, ebit ); }

template< size_t BIT, typename int_t >
inline void SetBit( int_t & v )
{ SetBits< BIT, BIT + 1 >( v ); }

template< typename int_t >
inline void SetBit( int_t & v, size_t bit )
{ v |= ((int_t)1<<bit); }

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline void ClearBits( int_t & v )
{ v &= SBitRangeTraits< int_t, SBIT, EBIT >::RMASK; }

template< typename int_t >
inline void ClearBitsBack( int_t & v, size_t bits )
{ v &= RMaskBack< int_t >( bits ); }

template< typename int_t >
inline void ClearBitsFront( int_t & v, size_t bits )
{ v &= RMaskFront< int_t >( bits ); }

template< typename int_t >
inline void ClearBits( int_t & v, size_t sbit, size_t ebit )
{ v &= RMask< int_t >( sbit, ebit ); }

template< size_t BIT, typename int_t >
inline void ClearBit( int_t & v )
{ ClearBits< BIT, BIT + 1 >( v ); }

template< typename int_t >
inline void ClearBit( int_t & v, size_t bit )
{ v &= (~((int_t)1<<bit)); }

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline void SelectBits( int_t & v )
{ v &= SBitRangeTraits< int_t, SBIT, EBIT >::MASK; }

template< typename int_t >
inline void SelectBitsBack( int_t & v, size_t bits )
{ v &= MaskBack< int_t >( bits ); }

template< typename int_t >
inline void SelectBitsFront( int_t & v, size_t bits )
{ v &= MaskFront< int_t >( bits ); }

template< typename int_t >
inline void SelectBits( int_t & v, size_t sbit, size_t ebit )
{ v &= Mask< int_t >( sbit, ebit ); }

template< size_t BIT, typename int_t >
inline void SelectBit( int_t & v )
{ SelectBits< BIT, BIT + 1 >( v ); }

template< typename int_t >
inline void SelectBit( int_t & v, size_t bit )
{ v &= ((int_t)1<<bit); }

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline int_t GetField( int_t v )
{ SelectBits< SBIT, EBIT >( v ); return (v>>SBIT); }

template< size_t BITS, typename int_t >
inline int_t GetFieldBack( int_t v )
{ SelectBits< 0, BITS >( v ); return v; }

template< typename int_t >
inline int_t GetFieldBack( size_t bits, int_t v )
{ SelectBitsBack( v, bits ); return v; }

template< typename int_t >
inline int_t GetField( size_t sbit, size_t ebit, int_t v )
{ v >>= sbit; SelectBitsBack( v, ebit - sbit ); return v; }

template< size_t BIT, typename int_t >
inline bool GetBit( int_t v )
{ return (bool)((v>>BIT)&0x1); }

template< typename int_t >
inline bool GetBit( size_t bit, int_t v )
{ return (bool)((v>>bit)&0x1); }

template< typename iter_t >
inline bool GetStreamBit( iter_t i, size_t bit )
{
    typedef typename SIterTraits< iter_t >::TValue TValue;
    static const size_t UBITS = BYTEBITS*sizeof( TValue );
    static const size_t SHIFT = SBinLog< UBITS >::VALUE;
    static const TValue MASK = SBitFieldTraits< TValue, SHIFT >::MASK;

    TValue unit( *(i + (bit>>SHIFT)) );
    bit = UBITS - (bit&MASK) - 1;
    return GetBit( bit, unit );
}

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline void AssignBits_Unsafe( int_t & d, int_t s )
{ ClearBits< SBIT, EBIT >( d ); d += s; }

template< typename int_t >
inline void AssignBits_Unsafe( size_t sbit, size_t ebit, int_t & d, int_t s )
{ ClearBits( d, sbit, ebit ); d += s; }

template< size_t SBIT, size_t EBIT, typename int_t >
inline void AssignBits( int_t & d, int_t s )
{ SelectBits< SBIT, EBIT >( s ); AssignBits_Unsafe< SBIT, EBIT >( d, s ); }

template< typename int_t >
inline void AssignBits( size_t sbit, size_t ebit, int_t & d, int_t s )
{ SelectBits( s, sbit, ebit ); AssignBits_Unsafe( sbit, ebit, d, s ); }

template< size_t BIT, typename int_t >
inline void AssignBit( int_t & v, bool bit )
{ if( bit ) SetBit< BIT >( v ); else ClearBit< BIT >( v ); }

template< typename int_t >
inline void AssignBit( size_t bit_n, int_t & v, bool bit )
{ if( bit ) SetBit( v, bit_n ); else ClearBit( v, bit_n ); }

//------------------------------------------------------------------------------
template< size_t SBIT, size_t EBIT, typename int_t >
inline void SetField_Unsafe( int_t & d, int_t s )
{ ClearBits< SBIT, EBIT >( d ); d += (s<<SBIT); }

template< typename int_t >
inline void SetField_Unsafe( size_t sbit, size_t ebit, int_t & d, int_t s )
{ ClearBits( d, sbit, ebit ); d += (s<<sbit); }

template< size_t BITS, typename int_t >
inline void SetFieldBack_Unsafe( int_t & d, int_t s )
{ ClearBits< 0, BITS >( d ); d += s; }

template< typename int_t >
inline void SetFieldBack_Unsafe( size_t bits, int_t & d, int_t s )
{ ClearBitsBack( d, bits ); d += s; }

template< size_t SBIT, size_t EBIT, typename int_t >
inline void SetField( int_t & d, int_t s )
{ SelectBits< 0, EBIT - SBIT >( s ); SetField_Unsafe< SBIT, EBIT >( d, s ); }

template< typename int_t >
inline void SetField( size_t sbit, size_t ebit, int_t & d, int_t s )
{ SelectBitsBack( s, ebit - sbit ); SetField_Unsafe( sbit, ebit, d, s ); }

template< size_t BITS, typename int_t >
inline void SetFieldBack( int_t & d, int_t s )
{ SelectBits< 0, BITS >( s ); SetFieldBack_Unsafe< BITS >( d, s ); }

template< typename int_t >
inline void SetFieldBack( size_t bits, int_t & d, int_t s )
{ SelectBitsBack( s, bits ); SetFieldBack_Unsafe( bits, d, s ); }

//------------------------------------------------------------------------------
template< size_t BITS, typename int_t >
inline void PushBits_Unsafe( int_t & d, int_t s )
{ d <<= BITS; SetFieldBack_Unsafe< BITS >( d, s ); }

template< typename int_t >
inline void PushBits_Unsafe( size_t bits, int_t & d, int_t s )
{ d <<= bits; SetFieldBack_Unsafe( bits, d, s ); }

template< size_t BITS, typename int_t >
inline void PushBits( int_t & d, int_t s )
{ SelectBits< 0, BITS >( s ); PushBits_Unsafe< BITS >( d, s ); }

template< typename int_t >
inline void PushBits( size_t bits, int_t & d, int_t s )
{ SelectBitsBack( s, 0 ); PushBits_Unsafe( bits, d, s ); }

//------------------------------------------------------------------------------
template< typename int_t >
inline int FirstSetBit_Left( int_t v )
{
    static const size_t WBITS = sizeof( int_t )*BYTEBITS;

    if( v == 0 ) return 1 + WBITS;
    size_t p( WBITS>>1 ), shift( p );

    while( p > 0 ) {
        p >>= 1;
        if( (int_t)(v>>shift) == 0 ) shift -= p; else shift += p;
    }

    return WBITS - shift - 1;
}

template< typename int_t >
inline int FirstSetBit_Right( int_t v )
{
    static const size_t WBITS = sizeof( int_t )*BYTEBITS;

    if( v == 0 ) return 1 + WBITS;
    size_t p( WBITS>>1 ), shift( p );

    while( p > 0 ) {
        p >>= 1;
        if( (int_t)(v<<shift) == 0 ) shift -= p; else shift += p;
    }

    return WBITS - shift - 1;
}

//------------------------------------------------------------------------------
template< typename int_t >
inline int CountBits( int_t v )
{
    int res( 0 );
    while( v != 0 ) { v &= v-1; ++res; }
    return res;
}

END_NS( common )
END_STD_SCOPES

#endif

