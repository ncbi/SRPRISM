/*  $Id: index_iterator.cpp 205423 2010-09-17 19:16:42Z morgulis $
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
 * File Description: high level code for .idx files traversal
 *
 */

#include <ncbi_pch.hpp>

#include "index_iterator.hpp"

#include <cassert>
#include <iostream>

#include "../common/util.hpp"
#include "../common/bits.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
const char * CIndexBase::IDX_MAP_SFX = ".map";
const char * CIndexBase::IDX_PROPER_SFX = ".idx";
const char * CIndexBase::IDX_REPMAP_SFX = ".rmp";

//------------------------------------------------------------------------------
CIndexIterator::CIndexIterator( const std::string & basename )
    : map_reader_( basename + IDX_MAP_SFX ),
      idx_reader_( basename + IDX_PROPER_SFX ),
      state_( NEW_PREFIX ), prefix_( 0 ), 
      init_( false ), special_( false ), end_( false ), start_( true )
{
}

//------------------------------------------------------------------------------
template< TStrand strand >
inline void CIndexIterator::ReadStrandDataSpecial(void)
{
    DBG_TRACE( "IDXITER: (special) strand: " << (int)strand );
    size_t npos;
    size_t next;
    TExtData & ext_data = ext_data_[strand];
    rv_pos_start_ = pos_.size();
    { Uint4 t; idx_reader_.NextWord( t ); next = (size_t)t; }
    { Uint4 t; idx_reader_.NextWord( t ); npos = (size_t)t; }
    ext_data.resize( next );
    pos_.reserve( rv_pos_start_ + npos );
    DBG_TRACE( "IDXITER: next: " << next << " npos: " << npos );

    for( size_t i = 0; i < next; ++i ) {
        SExtInfo & ext_info( ext_data[i] );
        { Uint4 t; idx_reader_.NextWord( t ); ext_info.ext   = (TUnit)t; }
        { Uint4 t; idx_reader_.NextWord( t ); ext_info.fnpos = (size_t)t; }
        { Uint4 t; idx_reader_.NextWord( t ); ext_info.rnpos = (size_t)t; }

        { 
            Uint4 t; idx_reader_.NextWord( t ); 
            ext_info.index = (size_t)t; 
        }

        DBG_TRACE( 
                "IDXITER: ext: " << std::hex << ext_info.ext << std::dec <<
                      " fnpos: " << ext_info.fnpos <<
                      " rnpos: " << ext_info.rnpos <<
                      " index: " << ext_info.index );
    }

    for( size_t i = 0; i < npos; ++i ) {
        Uint4 t; idx_reader_.NextWord( t );
        pos_.push_back( (TPos)t );
        DBG_TRACE( "IDXITER: pos: " << (TPos)t ); 
    }

    bytes_left_ -= sizeof( Uint4 )*(2 + npos + 4*next);
}

//------------------------------------------------------------------------------
void CIndexIterator::ReadDataSpecial(void)
{
    special_ = true; pos_.clear();
    ReadStrandDataSpecial< STRAND_FW >();
    ReadStrandDataSpecial< STRAND_RV >();
}

//------------------------------------------------------------------------------
void CIndexIterator::ReadData(void)
{
    special_ = false;
    rv_pos_start_ = 0;

    switch( f_bytes_ ) {
        case 0: break;

        case 1: 
            { 
                Uint1 t; idx_reader_.NextWord( t ); rv_pos_start_ = (size_t)t; 
                --bytes_left_; 
            }

            break;

        case 2:
            { 
                Uint2 t; idx_reader_.NextWord( t ); rv_pos_start_ = (size_t)t; 
                bytes_left_ -= 2; 
            }

            break;

        case 3:
            { 
                Uint4 t; idx_reader_.NextWord( t ); rv_pos_start_ = (size_t)t; 
                bytes_left_ -= 4; 
            }

            break;

        default: SRPRISM_ASSERT( false );
    }

    size_t rnpos = 0;

    switch( r_bytes_ ) {
        case 0: break;

        case 1: 
            { 
                Uint1 t; idx_reader_.NextWord( t ); rnpos = (size_t)t; 
                --bytes_left_; 
            }

            break;

        case 2:
            { 
                Uint2 t; idx_reader_.NextWord( t ); rnpos = (size_t)t; 
                bytes_left_ -= 2; 
            }

            break;

        case 3:
            { 
                Uint4 t; idx_reader_.NextWord( t ); rnpos = (size_t)t; 
                bytes_left_ -= 4; 
            }

            break;

        default: SRPRISM_ASSERT( false );
    }

    size_t npos = rv_pos_start_ + rnpos;
    pos_.resize( npos );
    DBG_TRACE( 
            "IDXITER: fwpos: " << rv_pos_start_ <<
                    " rvpos: " << rnpos <<
                     " npos: " << npos );

    for( size_t i = 0; i < npos; ++i ) {
        Uint4 t; idx_reader_.NextWord( t ); pos_[i] = (TPos)t;
        DBG_TRACE( "IDXITER: pos: " << pos_[i] );
    }

    bytes_left_ -= npos*sizeof( Uint4 );
}

//------------------------------------------------------------------------------
inline void CIndexIterator::CheckPalindrome(void)
{
    static const int BYTES = 1 + (PREFIX_LEN - 1)/COMPRESSION;
    static const TUnit MASK = common::SBitFieldTraits< TUnit, BYTEBITS >::MASK;

    TUnit rc_prefix = 0;

    for( int i = 0; i < BYTES; ++i ) {
        rc_prefix = (rc_prefix<<common::BYTEBITS) 
                  + seq::BYTE_RC[(prefix_>>(i*BYTEBITS))&MASK];
    }

    palindrome_ = (prefix_ == rc_prefix);
}

//------------------------------------------------------------------------------
bool CIndexIterator::Next(void)
{
    if( end_ ) { state_ = END_OF_INDEX; return false; }

    switch( state_ ) {
        case NEW_PREFIX:
            {
                special_ = false;

                { 
                    Uint4 t; 
                    idx_reader_.NextWord( t ); 
                    bytes_left_ = (size_t)t; 
                    DBG_TRACE( "IDXITER: bytes left: " << bytes_left_ );
                }
                
                if( idx_reader_.Eof() ) state_ = END_OF_INDEX;
                else { state_ = IN_PREFIX; return Next(); }

                return false;
            }
        
        case IN_PREFIX:
            {
                prefix_ = Map2Unit( map_reader_.Unit() );
                TUnit sfx; 
                { Uint2 t; idx_reader_.NextWord( t ); sfx = (TUnit)t; }
                bytes_left_ -= EXACT_BYTES;
                r_bytes_ = (size_t)GetField( 
					R_BYTES_START_BIT, R_BYTES_END_BIT, sfx );
                f_bytes_ = (size_t)GetField( 
					F_BYTES_START_BIT, F_BYTES_END_BIT, sfx );
                sfx >>= SFX_START_BIT;
                prefix_ += sfx;
                DBG_TRACE( 
                        "IDXITER r_bytes: " << r_bytes_ <<
                               " f_bytes: " << f_bytes_ <<
                                " prefix: " << std::hex << prefix_ << 
                         std::dec );
                CheckPalindrome();
                if( f_bytes_ == 0 && r_bytes_ == 0 ) ReadDataSpecial();
                else ReadData();

                if( bytes_left_ == 0 ) {
                    if( map_reader_.Advance() ) {
                        state_ = NEW_PREFIX;
                        idx_reader_.FF( map_reader_.Offset() );
                    }
                    else {
                        end_ = true;
                        return false;
                    }
                }

                return true;
            }

        default: SRPRISM_ASSERT( false );

        return true;
    }
}

//------------------------------------------------------------------------------
bool CIndexIterator::Seek( TUnit prefix )
{
    if( !start_ && prefix <= prefix_ ) return true;
    start_ = false;
    if( !map_reader_.Seek( Unit2Map( prefix ) ) ) return false;
    idx_reader_.FF( map_reader_.Offset() );
    if( Unit2Map( prefix_ ) < map_reader_.Unit() ) state_ = NEW_PREFIX;
    if( !init_ ) { Next(); init_ = true; }

    while( state_ != END_OF_INDEX && prefix_ < prefix ) {
        if( !Next() ) break;
    }

    return (state_ != END_OF_INDEX);
}

END_NS( srprism )
END_STD_SCOPES

