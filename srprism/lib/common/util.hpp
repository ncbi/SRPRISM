/*  $Id: util.hpp 205411 2010-09-17 17:42:11Z morgulis $
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
 * File Description: miscelanious utilities
 *
 */

#ifndef __AM_COMMON_UTIL_HPP__
#define __AM_COMMON_UTIL_HPP__

#include "binfile.hpp"

START_STD_SCOPES
START_NS( common )

//
// function objects
//
template< typename lhs_t, typename rhs_t >
struct CPlus
{
    CPlus( const rhs_t & c = rhs_t() ) : c_( c ) {}

    const lhs_t operator()( const lhs_t & lhs, const rhs_t & rhs ) const
    { return lhs + rhs + c_; }

    private: const rhs_t & c_;
};

//
// algorithms
//
template< typename func_t, typename val_t, typename iter_t >
const val_t rfold( const func_t & f, const val_t & v, iter_t s, iter_t e )
{
    val_t res( v );
    for( iter_t i = s; i != e; ++i ) res = f( res, *i );
    return res;
}

//
// host endianness
//
struct SEndianness
{
    typedef Uint1 TTYPE;
    static const TTYPE LITTLE;
    static const TTYPE BIG;
};

extern SEndianness::TTYPE Endianness();

//
// reading integers from binary files
//
template< typename int_t >
bool ReadIntDirect( int_t & res, CReadBinFile & file )
{ return file.Read( (char *)&res, sizeof( int_t ), true ) == sizeof( int_t ); }

template< typename int_t >
bool ReadInt( int_t & res, CReadBinFile & file, bool endian_match )
{
    if( endian_match ) return ReadIntDirect( res, file );

    res = 0;
    char buf[sizeof(int_t)];

    if( file.Read( buf, sizeof( int_t ), true ) != sizeof( int_t ) ) {
        return false;
    }

    for( Uint1 i = 0; i < sizeof( int ); ++i ) res += (buf[i]<<(i<<3));
    return true;
}

END_NS( common )
END_STD_SCOPES

#endif
