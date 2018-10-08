/*  $Id: seqdef.cpp 205414 2010-09-17 17:59:42Z morgulis $
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

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <cassert>
#include <algorithm>

#include "../common/trace.hpp"
#include "seqdef.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

const size_t BYTE_RC_LEN 
    = common::SBitFieldTraits< size_t, common::BYTEBITS >::MAX + 1;
common::Uint1 BYTE_RC[BYTE_RC_LEN];

TAlphabetString SCodingTraits_Base< CODING_IUPACNA >::ALPHABET_STRING;
TAlphabetArray SCodingTraits_Base< CODING_NCBI2NA >::ALPHABET_ARRAY =
{'\000','\001','\002','\003'};
TAlphabetArray SCodingTraits_Base< CODING_NCBI4NA >::ALPHABET_ARRAY =
{'\000','\001','\002','\003','\004','\005','\006','\007',
 '\010','\011','\012','\013','\014','\015','\016','\017'};

TAlphabetString SCodingTraits_Base< CODING_IUPACNA >::RC_ALPHABET_STRING;
TAlphabetArray SCodingTraits_Base< CODING_NCBI2NA >::RC_ALPHABET_ARRAY =
{'\003','\002','\001','\000'};
TAlphabetArray SCodingTraits_Base< CODING_NCBI4NA >::RC_ALPHABET_ARRAY =
{'\000','\010','\004','\014','\002','\012','\006','\016',
 '\001','\011','\005','\015','\003','\013','\007','\017'};

TAlphabetString SCodingTraits_Base< CODING_IUPACNA >::AMBIG_STRING;
TAlphabetArray SCodingTraits_Base< CODING_NCBI2NA >::AMBIG_ARRAY =
{ 'f', 'f', 'f', 'f' };
TAlphabetArray SCodingTraits_Base< CODING_NCBI4NA >::AMBIG_ARRAY =
{ 't', 'f', 'f', 't', 'f', 't', 't', 't', 
  'f', 't', 't', 't', 't', 't', 't', 't' };

template< TCoding coding > TAlphabet SCodingTraits< coding >::ALPHABET;

template< TCoding coding > bool 
SCodingTraits< coding >::NAMBIG[SCodingTraits< coding >::NUM_VALUES];

template< TCoding coding > TLetter 
SCodingTraits< coding >::RC[SCodingTraits< coding >::NUM_VALUES];

template<> template<>
TLetter SCodingTraits< CODING_NCBI4NA >::Recode< CODING_NCBI4NA >(
	TLetter l )
{ return l; }

template<> template<>
TLetter SCodingTraits< CODING_IUPACNA >::Recode< CODING_NCBI2NA >(
        TLetter l )
{
    switch( l ) {
        case 'A': case 'a': return (TLetter)0;
        case 'C': case 'c': return (TLetter)1;
        case 'G': case 'g': return (TLetter)2;
        case 'T': case 't': return (TLetter)3;
    }

    return 0;
}

template<> template<>
TLetter SCodingTraits< CODING_IUPACNA >::Recode< CODING_NCBI4NA >(
        TLetter l )
{
    switch( l ) {
        case '-': return (TLetter)0;
        case 'A': case 'a': return (TLetter)1;
        case 'C': case 'c': return (TLetter)2;
        case 'M': case 'm': return (TLetter)3;
        case 'G': case 'g': return (TLetter)4;
        case 'R': case 'r': return (TLetter)5;
        case 'S': case 's': return (TLetter)6;
        case 'V': case 'v': return (TLetter)7;
        case 'T': case 't': return (TLetter)8;
        case 'W': case 'w': return (TLetter)9;
        case 'Y': case 'y': return (TLetter)10;
        case 'H': case 'h': return (TLetter)11;
        case 'K': case 'k': return (TLetter)12;
        case 'D': case 'd': return (TLetter)13;
        case 'B': case 'b': return (TLetter)14;
        case 'N': case 'n': return (TLetter)15;
    }

    return 0;
}

template<> template<>
TLetter SCodingTraits< CODING_NCBI2NA >::Recode< CODING_IUPACNA >(
        TLetter l )
{
    switch( l ) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
    }

    SRPRISM_ASSERT( false );
    return 0;
}

template<> template<>
TLetter SCodingTraits< CODING_NCBI2NA >::Recode< CODING_NCBI4NA > (
        TLetter l )
{
    switch( l ) {
        case 0: return 1;
        case 1: return 2;
        case 2: return 4;
        case 3: return 8;
    }

    SRPRISM_ASSERT( false );
    return 0;
}

template<> template<>
TLetter SCodingTraits< CODING_NCBI2NA >::Recode< CODING_NCBI2NA > (
        TLetter l )
{ return l; }

template<> template<>
TLetter SCodingTraits< CODING_NCBI4NA >::Recode< CODING_IUPACNA >(
        TLetter l )
{
    switch( l ) {
        case 0:  return '-';
        case 1:  return 'A';
        case 2:  return 'C';
        case 3:  return 'M';
        case 4:  return 'G';
        case 5:  return 'R';
        case 6:  return 'S';
        case 7:  return 'V';
        case 8:  return 'T';
        case 9:  return 'W';
        case 10: return 'Y';
        case 11: return 'H';
        case 12: return 'K';
        case 13: return 'D';
        case 14: return 'B';
        case 15: return 'N';
    }

    SRPRISM_ASSERT( false );
    return 0;
}

template<> template<>
TLetter SCodingTraits< CODING_NCBI4NA >::Recode< CODING_NCBI2NA >(
        TLetter l )
{
    switch( l ) {
        case 1: return 0;
        case 2: return 1;
        case 4: return 2;
        case 8: return 3;
        default: return 0;
    }

    SRPRISM_ASSERT( false );
}

template< TCoding coding > void InitForCoding(void)
{
    DBG_TRACE( "InitForCoding< " << (int)coding << ">" );
    typedef SCodingTraits_Base< coding > TBase;
    typedef SCodingTraits< coding > TType;

    for( size_t i = 0; i < TType::NUM_VALUES; ++i ) {
        TType::ALPHABET.push_back( TBase::ALPHABET_ARRAY[i] );
    }

    for( TAlphabet::const_iterator i = TType::ALPHABET.begin();
            i != TType::ALPHABET.end(); ++i ) {
        if( TBase::AMBIG_ARRAY[i - TType::ALPHABET.begin()] == 'f' ) {
            (TType::NAMBIG)[*i] = true;
        }

        (TType::RC)[*i] = 
            TBase::RC_ALPHABET_ARRAY[i - TType::ALPHABET.begin()];
    }
}

template<> void InitForCoding< CODING_IUPACNA >(void)
{
    static const TCoding coding = CODING_IUPACNA;
    DBG_TRACE( "InitForCoding< " << (int)coding << ">" );
    typedef SCodingTraits_Base< coding > TBase;
    typedef SCodingTraits< coding > TType;

    std::copy( 
        TBase::ALPHABET_STRING.begin(), TBase::ALPHABET_STRING.end(),
        std::back_inserter( TType::ALPHABET ) );

    for( TAlphabet::const_iterator i = TType::ALPHABET.begin();
            i != TType::ALPHABET.end(); ++i ) {
        if( TBase::AMBIG_STRING[i - TType::ALPHABET.begin()] == 'f' ) {
            (TType::NAMBIG)[*i] = true;
        }

        (TType::RC)[*i] = 
            TBase::RC_ALPHABET_STRING[i - TType::ALPHABET.begin()];
    }
}

void InitCoding(void)
{
    SCodingTraits_Base< CODING_IUPACNA >::ALPHABET_STRING = 
        "ACGTURYKMSWBDHVN-acgturykmswbdhvn";
    SCodingTraits_Base< CODING_IUPACNA >::AMBIG_STRING = 
        "fffftttttttttttttfffftttttttttttt";
    SCodingTraits_Base< CODING_IUPACNA >::RC_ALPHABET_STRING = 
        "TGCAAYRMKSWVHDBN-tgcaayrmkswvhdbn";

    InitForCoding< CODING_IUPACNA >();
    InitForCoding< CODING_NCBI2NA >();
    InitForCoding< CODING_NCBI4NA >();

    {
        static const TCoding CODING = CODING_NCBI2NA;
        static const TSeqSize NUM_LETTERS = 
            SCodingTraits< CODING >::PACK_FACTOR;
        typedef common::Uint1 TUnit;

        for( size_t i = 0; i < BYTE_RC_LEN; ++i ) {
            for( TSeqSize j = 0; j < NUM_LETTERS; ++j ) {
                TLetter l = GetLetter< CODING >( (TUnit)i, j );
                SetLetter< CODING >( 
                        BYTE_RC[i], NUM_LETTERS - j - 1, 
                        SCodingTraits< CODING >::RC[l] );
            }
        }
    }
}

const std::string CS_ALPHABET_STRING = "0123";

static const TLetter CS_MATRIX[4][4] = {
    { 0, 1, 2, 3 },
    { 1, 0, 3, 2 },
    { 2, 3, 0, 1 },
    { 3, 2, 1, 0 }
};

std::string Color2IUPACNA( const std::string & c_str, TLetter lead )
{
    std::string res;

    for( std::string::const_iterator ci( c_str.begin() ); 
            ci != c_str.end(); ++ci ) {
        int offset( 0 );

        switch( *ci ) {
            case '0': offset = 0; break;
            case '1': offset = 1; break;
            case '2': offset = 2; break;
            case '3': offset = 3; break;
            default: SRPRISM_ASSERT( false ); break;
        }

        lead = CS_MATRIX[lead][offset];
        res.push_back( SCodingTraits< CODING_NCBI2NA >::Recode< 
                CODING_IUPACNA >( lead ) );
    }

    return res;
}

std::string IUPACNA2Color( const std::string & str, TLetter lead )
{
    std::string res; res.push_back( lead );
    lead = SCodingTraits< CODING_IUPACNA >::Recode< CODING_NCBI2NA >( lead );

    for( std::string::const_iterator ci( str.begin() );
            ci != str.end(); ++ci ) {
        TLetter l( SCodingTraits< CODING_IUPACNA >::Recode< 
                CODING_NCBI2NA >( *ci ) );
        TLetter cl( CS_MATRIX[lead][l] );
        lead = l;

        switch( cl ) {
            case 0: res.push_back( '0' ); break;
            case 1: res.push_back( '1' ); break;
            case 2: res.push_back( '2' ); break;
            case 3: res.push_back( '3' ); break;
            default: SRPRISM_ASSERT( false ); break;
        }
    }

    return res;
}

END_NS( seq )
END_STD_SCOPES
