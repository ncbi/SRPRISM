/*  $Id: query_data.cpp 348993 2012-01-06 14:51:45Z morgulis $
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
 * File Description: unique'ed query data structure
 *
 */

#include <ncbi_pch.hpp>

#include "../seq/seqdef.hpp"
#include "query_data.hpp"
#include "seqiter.hpp"

START_STD_SCOPES
START_NS( srprism )

USE_NS( seq )
USE_NS( common )

static const int DIV_SHIFT = 2;

//------------------------------------------------------------------------------
bool CQueryData::ComputeHC_NormalFixedLen( 
        const CRMap & rmap, TSeqSize sa_start, TSeqSize sa_end,
        TSeqSize start, TSeqSize end, 
        common::Uint8 min_badness,
        const Uint8 * bvec, const Uint8 * mvec, const Uint8 * avec,
        TWord & hc, TWord & mw, Uint8 & badness )
{
    if( !IsShort() && 
            start - sa_start < HASH_LEN &&
            sa_end - start < HASH_LEN + HASH_LEN ) {
        return false;
    }

    for( size_t i( start ); i < end; i += HASH_LEN ) {
        badness += bvec[i - sa_start];
    }

    if( badness >= min_badness ) return false;

    int n_hashes( (end - start)/HASH_LEN );
    TSeqSize start_off( start - sa_start ), end_off( end - sa_start );

    for( size_t i( start_off ); i < end_off; i += HASH_LEN ) {
        if( i < HASH_LEN ) {
            if( sa_end - sa_start - i < HASH_LEN + HASH_LEN ) {
                common::AssignBit( 
                        SA_EXT_DIR_MASK_START_BIT + (i - start_off)/HASH_LEN, 
                        mw, (i < sa_end - sa_start - i - HASH_LEN) );
            }
            else {
                common::SetBit( 
                        mw, 
                        SA_EXT_DIR_MASK_START_BIT + (i - start_off)/HASH_LEN );
            }
        }
        else if( sa_end - sa_start - i >= HASH_LEN + HASH_LEN ) {
            common::AssignBit( 
                    SA_EXT_DIR_MASK_START_BIT + (i - start_off)/HASH_LEN, mw, 
                    (mvec[i] == 1) );
        }
    }

    int n_ambigs( 0 );

    if( IsAmbig() ) {
        for( size_t i( start ); i < end; i += HASH_LEN ) {
            n_ambigs += avec[i - sa_start];
        }
    }
    
    if( n_hashes <= n_ambigs ) return false;
    common::SetField< SA_OFF_START_BIT, SA_OFF_END_BIT >( hc, (TWord)start );
    common::SetField< SA_N_HASHES_START_BIT, SA_N_HASHES_END_BIT >(
            hc, (TWord)(n_hashes - 1) );
    common::SetBit< SA_LONG_BIT >( hc );
    common::AssignBit< SA_AMBIG_BIT >( hc, (n_ambigs > 0) );
    return true;
}

//------------------------------------------------------------------------------
bool CQueryData::ComputeHC_Normal( 
        const CRMap & rmap, TSeqSize sa_start, TSeqSize sa_end, int n_err )
{
    TSeqSize sa_len( (n_err + 1)*HASH_LEN );
    TWord hc( 0 ), mw( 0 );
    Uint8 min_badness( SIntTraits< Uint8 >::MAX );
    TSeqSize end( sa_end - sa_start - HASH_LEN + 1 );
    std::vector< Uint8 > bvec( end, 0 ), mvec( end, 0 ), avec( end, 0 );
    CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si(
            Data(), sa_start, sa_end );

    for( TSeqSize i( 0 ); i < end; ++i, si += 1 ) {
        bvec[i] = GetWordBadness( rmap, (*si).first );
    }

    for( TSeqSize i( HASH_LEN ); i + HASH_LEN < end; ++i ) {
        mvec[i] = (bvec[i - HASH_LEN] <= bvec[i + HASH_LEN]) ? 0 : 1;
    }

    if( IsAmbig() ) {
        for( TSeqSize i( 0 ); i < end; ++i ) {
            avec[i] = NHashAmbigs( sa_start + i );
        }
    }

    end = sa_end - sa_len + 1;

    for( TSeqSize i( sa_start ); i < end; ++i ) {
        TWord hc_1( 0 ), mw_1( 0 );
        Uint8 b( 0 );

        if( ComputeHC_NormalFixedLen( 
                    rmap, sa_start, sa_end, i, i + sa_len, min_badness, 
                    &bvec[0], &mvec[0], &avec[0], 
                    hc_1, mw_1, b ) ) {
            if( b < min_badness ) {
                min_badness = b;
                hc = hc_1;
                mw = mw_1;
            }
        }
    }

    if( min_badness != SIntTraits< Uint8 >::MAX ) {
        *raw_data_ = hc;
        *(raw_data_ + 1) = mw;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------
bool CQueryData::ComputeHC_Short( 
        const CRMap & rmap, TSeqSize sa_start, TSeqSize sa_end )
{
    int n_ambigs( 0 );
    TWord & hc( *raw_data_ ), & mw( *(raw_data_ + 1) );
    TSeqSize s1( sa_end - HASH_LEN );

    if( IsAmbig() ) {
        n_ambigs = NAmbigs( sa_start, sa_end );
        if( n_ambigs > 1 ) return false;
        AssignBit< SA_AMBIG_BIT >( hc, (n_ambigs > 0) );
    }

    Uint8 bad_0, bad_1( 0 );
    CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si( 
            Data(), sa_start, sa_end );
    TWord h0( (*si).first ), h1( 0 );
    bad_0 = GetWordBadness( rmap, h0 );
    SetField< SA_OFF_START_BIT, SA_OFF_END_BIT >( hc, sa_start );
    SetField< SA_N_HASHES_START_BIT, SA_N_HASHES_END_BIT >( hc, (TWord)1 );
    SetField< SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >( 
            hc, s1 - sa_start );
    SetBit< SA_EXT_DIR_MASK_START_BIT >( mw );

    if( s1 > sa_start ) {
        si += (s1 - sa_start);
        h1 = (*si).first;
        bad_1 = GetWordBadness( rmap, h1 );

        if( n_ambigs > 0 && !(HashAmbig( 0 ) && HashAmbig( 1 )) ) {
            SetField< 
                SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >(
                        hc, (TWord)3 );
        }
        else {
            TWord bu_hash_idx( bad_0 <= bad_1 ? 0 : 1 );
            SetField<
                SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >(
                        hc, bu_hash_idx );
            AssignBit< SA_BLOWUP_DIR_BIT >( hc, (bad_0 <= bad_1) );
        }
    }
    else {
        SetField< SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >(
                hc, (TWord)0 );
        SetBit< SA_BLOWUP_DIR_BIT >( hc );
    }

    return true;
}

//------------------------------------------------------------------------------
bool CQueryData::ComputeHC_Med( 
        const CRMap & rmap, TSeqSize sa_start, TSeqSize sa_end )
{
    int n_ambigs( 0 );
    TWord & hc( *raw_data_ ), & mw( *(raw_data_ + 1) );

    if( IsAmbig() ) {
        n_ambigs = NAmbigs( sa_start, sa_end );
        if( n_ambigs > 2 ) return false;
        AssignBit< SA_AMBIG_BIT >( hc, (n_ambigs > 0) );
    }

    TSeqSize s0( sa_start ), s1( sa_end - HASH_LEN - HASH_LEN ),
             s2( sa_start + HASH_LEN ), s3( s1 + HASH_LEN );
    Uint8 b0, b1, b2, b3;
    CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si( 
            Data(), sa_start, sa_end );
    b0 = GetWordBadness( rmap, (*si).first ); si += s1 - s0;
    b1 = GetWordBadness( rmap, (*si).first ); si += s2 - s1;
    b2 = GetWordBadness( rmap, (*si).first ); si += s3 - s2;
    b3 = GetWordBadness( rmap, (*si).first );
    SetField< SA_OFF_START_BIT, SA_OFF_END_BIT >( hc, sa_start );
    SetField< SA_N_HASHES_START_BIT, SA_N_HASHES_END_BIT >( hc, (TWord)2 );
    SetField< SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >( 
            hc, s3 - sa_start - HASH_LEN );
    SetBit< SA_EXT_DIR_MASK_START_BIT >( mw );
    AssignBit< SA_EXT_DIR_MASK_START_BIT + 1 >( mw, (b1 <= b2) );

    {
        int idx; bool right_bu_ext;

        if( b1 <= b2 ) {
            if( s1 == s0 || b0 <= b1 ) { idx = 0; right_bu_ext = true; }
            else{ idx = 1; right_bu_ext = false; }
        }
        else {
            if( s2 == s3 || b3 > b2 ) { idx = 2; right_bu_ext = false; }
            else { idx = 1; right_bu_ext = true; }
        }

        SetField< SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                hc, (TWord)idx );
        AssignBit< SA_BLOWUP_DIR_BIT >( hc, right_bu_ext );
    }

    if( n_ambigs > 0 ) {
        int h0_ambigs( NHashAmbigs( s0 ) ),
            h2_ambigs( NHashAmbigs( s3 ) ),
            h1_ambigs( (b1 <= b2) ? NHashAmbigs( s1 ) : NHashAmbigs( s2 ) );

        if( h0_ambigs > 1 || h1_ambigs > 1 || h2_ambigs > 1 ) {
            SetField< 
                SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >(
                        hc, (TWord)3 );
        }
        else if( n_ambigs == 2 ) {
            if( h0_ambigs == 0 || h1_ambigs == 0 || h2_ambigs == 0 ) {
                SetField< 
                    SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                            hc, (TWord)3 );
            }
        }
        else if( b1 > b2 ) {
            if( h0_ambigs > 0 ) {
                SetField< 
                    SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                            hc, (TWord)0 );
                SetBit< SA_BLOWUP_DIR_BIT >( hc );
            }
            else if( h1_ambigs == 0 || h2_ambigs == 0 ) {
                SetField< 
                    SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                            hc, (TWord)3 );
            }
        }
        else {
            if( h2_ambigs > 0 ) {
                SetField< 
                    SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                            hc, (TWord)2 );
                ClearBit< SA_BLOWUP_DIR_BIT >( hc );
            }
            else if( h1_ambigs == 0 || h0_ambigs == 0 ) {
                SetField< 
                    SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >( 
                            hc, (TWord)3 );
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
CQueryData::CBlowUpIterator_Base::CBlowUpIterator_Base( 
        const CQueryData & q, TSeqSize hash_start, 
        TSeqSize start, TSeqSize end, bool right_ext )
    : q_( &q ), hash_start_( hash_start ), start_( start ), end_( end ),
      type_1_( SErrType::E ), type_2_( SErrType::E ),
      right_ext_( right_ext ),
      done_( false ), two_err_( false ),
      curr_pos_adj_( right_ext_ ? 1 : 0 )
{
    CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si( 
            q_->Data(), hash_start_, hash_start_ + HASH_LEN );
    orig_hash_ = (*si).first;
}

//------------------------------------------------------------------------------
const seq::TAlphabet & CQueryData::CBlowUpIterator_Base::alphabet = 
    SCodingTraits< CQueryData::CBlowUpIterator_Base::CODING >::ALPHABET;

//------------------------------------------------------------------------------
CQueryData::CBlowUpIterator_M::CBlowUpIterator_M( 
        const CQueryData & q, TSeqSize hash_start, 
        TSeqSize start, TSeqSize end, bool right_ext, bool ambig )
    : CBlowUpIterator_Base( q, hash_start, start, end, right_ext ), 
      alpha_( alphabet.begin() ), ambig_( ambig )
{
    if( start >= end ) { done_ = true; return; }
    type_1_ = SErrType::M;
    curr_pos_adj_ = 1;
    curr_pos_ = start_;
    qletter_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
    if( *alpha_ == qletter_ && !ambig_ ) ++alpha_;
    hash_ = orig_hash_;
    SetLetter< CODING >( hash_, curr_pos_, *alpha_ );
}

//------------------------------------------------------------------------------
bool CQueryData::CBlowUpIterator_M::Next( void )
{
    if( End() ) return false;

    do { ++alpha_; } 
    while( alpha_ != alphabet.end() && *alpha_ == qletter_ && !ambig_ );

    if( alpha_ == alphabet.end() ) {
        if( ++curr_pos_ >= end_ ) { done_ = true; return false; }
        qletter_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
        alpha_ = alphabet.begin();
        if( *alpha_ == qletter_ && !ambig_ ) ++alpha_;
        hash_ = orig_hash_;
    }

    SetLetter< CODING >( hash_, curr_pos_, *alpha_ ); 
    return true; 
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_I::RemoveAndShiftLeft( void )
{
    hash_ = orig_hash_;
    TWord t( SelectLetters< CODING >( curr_pos_ + 1, HASH_LEN, hash_ ) );
    t <<= LETTER_BITS;
    AssignLetters_Unsafe< CODING >( curr_pos_, HASH_LEN - 1, hash_, t );
    TLetter ll( GetStreamLetter< CODING >( q_->Data(), hash_start_ + HASH_LEN ) );
    SetLetter< CODING >( hash_, HASH_LEN - 1, ll );
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_I::RemoveAndShiftRight( void )
{
    hash_ = orig_hash_;
    TWord t( SelectLetters< CODING >( 0, curr_pos_ - 1, hash_ ) );
    t >>= LETTER_BITS;
    AssignLetters_Unsafe< CODING >( 1, curr_pos_, hash_, t );
    TLetter ll( GetStreamLetter< CODING >( 
                q_->Data(), hash_start_ - 1 ) );
    SetLetter< CODING >( hash_, 0, ll );
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_I::SearchRight( void )
{
    prev_letter_ = qletter_1_;
    qletter_1_ = GetLetter< CODING >( orig_hash_, curr_pos_ );

    while( curr_pos_ < end_ && prev_letter_ == qletter_1_ ) {
        ++curr_pos_;
        prev_letter_ = qletter_1_;
        qletter_1_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
    }
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_I::SearchLeft( void )
{
    prev_letter_ = qletter_1_;
    qletter_1_ = GetLetter< CODING >( orig_hash_, curr_pos_ - 1 );

    while( curr_pos_ > start_ && prev_letter_ == qletter_1_ ) {
        --curr_pos_;
        prev_letter_ = qletter_1_;
        qletter_1_ = GetLetter< CODING >( orig_hash_, curr_pos_ - 1 );
    }
}

//------------------------------------------------------------------------------
CQueryData::CBlowUpIterator_I::CBlowUpIterator_I( 
        const CQueryData & q, TSeqSize hash_start, 
        TSeqSize start, TSeqSize end, bool right_ext, bool ambig )
    : CBlowUpIterator_Base( q, hash_start, start, end, right_ext ),
      qletter_1_( 0xFF ), qletter_2_( 0xFF )
{
    if( q.Len() <= HASH_LEN ) { done_ = true; return; }
    if( start >= end ) { done_ = true; return; }
    hash_ = orig_hash_;
    type_1_ = SErrType::I;

    if( right_ext ) {
        curr_pos_ = start_;

        if( !ambig && hash_start_ > 0 ) {
            qletter_1_ = GetStreamLetter< CODING >( 
                    q_->Data(), hash_start_ + start_ - 1 );
            SearchRight();
        }
        else {
            qletter_1_ = GetStreamLetter< CODING >(
                    q_->Data(), hash_start_ + start_ );
        }

        if( curr_pos_ == end_ ) { done_ = true; return; }
        RemoveAndShiftLeft();
    }
    else {
        curr_pos_ = end_;

        if( !ambig && hash_start_ + end_ < q_->Len() ) {
            qletter_1_ = GetStreamLetter< CODING >(
                    q_->Data(), hash_start_ + end_ );
            SearchLeft();
        }
        else {
            qletter_1_ = GetStreamLetter< CODING >(
                    q_->Data(), hash_start_ + end_ - 1 );
        }

        if( curr_pos_ == start_ ) { done_ = true; return; }
        RemoveAndShiftRight();
    }
}

//------------------------------------------------------------------------------
bool CQueryData::CBlowUpIterator_I::Next( void )
{
    if( End() ) return false;

    if( right_ext_ ) {
        if( !two_err_ ) {
            if( q_->IsShort() || end_ - start_ < HASH_LEN ) {
                if( ++curr_pos_ == end_ ) { done_ = true; return false; }
                SearchRight();
                if( curr_pos_ == end_ ) { done_ = true; return false; }
                RemoveAndShiftLeft();
            }
            else {
                two_err_ = true;
                alpha_ = alphabet.begin();
                type_2_ = SErrType::M;
                qletter_2_ = GetStreamLetter< CODING >( 
                        q_->Data(), hash_start_ + end_ );

                if( *alpha_ == qletter_2_ ) {
                    if( hash_start_ + end_ + 1 < q_->Len() ) {
                        TLetter ll( GetStreamLetter< CODING >(
                                q_->Data(), hash_start_ + end_ + 1 ) );
                        SetLetter< CODING >( hash_, end_ - 1, ll );
                        type_2_ = SErrType::I;
                    }
                    else ++alpha_;
                }

                if( type_2_ == SErrType::M ) {
                    SetLetter< CODING >( hash_, end_ - 1, *alpha_ );
                }
            }
        }
        else {
            if( ++alpha_ == alphabet.end() ) {
                if( ++curr_pos_ == end_ ) { done_ = true; return false; }
                two_err_ = false;
                type_2_ = SErrType::E;
                SearchRight();
                if( curr_pos_ == end_ ) { done_ = true; return false; }
                RemoveAndShiftLeft();
            }
            else {
                type_2_ = SErrType::M;
                qletter_2_ = GetStreamLetter< CODING >( 
                        q_->Data(), hash_start_ + end_ );

                if( *alpha_ == qletter_2_ ) {
                    if( hash_start_ + end_ + 1 < q_->Len() ) {
                        TLetter ll( GetStreamLetter< CODING >(
                                q_->Data(), hash_start_ + end_ + 1 ) );
                        SetLetter< CODING >( hash_, end_ - 1, ll );
                        type_2_ = SErrType::I;
                    }
                    else ++alpha_;
                }

                if( type_2_ == SErrType::M ) {
                    SetLetter< CODING >( hash_, end_ - 1, *alpha_ );
                }
            }
        }
    }
    else {
        if( !two_err_ ) {
            if( q_->IsShort() || end_ - start_ < HASH_LEN ) {
                if( --curr_pos_ == start_ ) { done_ = true; return false; }
                SearchLeft();
                if( curr_pos_ == start_ ) { done_ = true; return false; }
                RemoveAndShiftRight();
            }
            else {
                two_err_ = true;
                alpha_ = alphabet.begin();
                type_1_ = SErrType::M;
                type_2_ = SErrType::I;
                qletter_2_ = qletter_1_;
                qletter_1_ = GetStreamLetter< CODING >( 
                        q_->Data(), hash_start_ + start_ - 1 );

                if( *alpha_ == qletter_1_ ) {
                    if( hash_start_ + start_ > 1 ) {
                        TLetter ll( GetStreamLetter< CODING >(
                                q_->Data(), hash_start_ + start_ - 2 ) );
                        SetLetter< CODING >( hash_, start_, ll );
                        type_1_ = SErrType::I;
                    }
                    else ++alpha_;
                }

                if( type_1_ == SErrType::M ) {
                    SetLetter< CODING >( hash_, start_, *alpha_ );
                }
            }
        }
        else {
            if( ++alpha_ == alphabet.end() ) {
                qletter_1_ = qletter_2_;
                type_1_ = SErrType::I;
                type_2_ = SErrType::E;
                two_err_ = false;
                if( --curr_pos_ == start_ ) { done_ = true; return false; }
                SearchLeft();
                if( curr_pos_ == start_ ) { done_ = true; return false; }
                RemoveAndShiftRight();
            }
            else {
                type_1_ = SErrType::M;
                qletter_1_ = GetStreamLetter< CODING >( 
                        q_->Data(), hash_start_ + start_ - 1 );

                if( *alpha_ == qletter_1_ ) {
                    if( hash_start_ + start_ > 1 ) {
                        TLetter ll( GetStreamLetter< CODING >(
                                q_->Data(), hash_start_ + start_ - 2 ) );
                        SetLetter< CODING >( hash_, start_, ll );
                        type_1_ = SErrType::I;
                    }
                    else ++alpha_;
                }

                if( type_1_ == SErrType::M ) {
                    SetLetter< CODING >( hash_, start_, *alpha_ );
                }
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_D::InsertAndShiftRight( void )
{
    TWord t( SelectLetters< CODING >( curr_pos_, end_ - 1, orig_hash_ ) );
    t >>= LETTER_BITS;
    AssignLetters< CODING >( curr_pos_ + 1, end_, hash_, t );
    SetLetter< CODING >( hash_, curr_pos_, *alpha_ );
}

//------------------------------------------------------------------------------
void CQueryData::CBlowUpIterator_D::InsertAndShiftLeft( void )
{
    TWord t( SelectLetters< CODING >( start_ + 1, curr_pos_, orig_hash_ ) );
    t <<= LETTER_BITS;
    AssignLetters< CODING >( start_, curr_pos_ - 1, hash_, t );
    SetLetter< CODING >( hash_, curr_pos_ - 1, *alpha_ );
}

//------------------------------------------------------------------------------
CQueryData::CBlowUpIterator_D::CBlowUpIterator_D( 
        const CQueryData & q, TSeqSize hash_start, 
        TSeqSize start, TSeqSize end, bool right_ext, bool ambig )
    : CBlowUpIterator_Base( q, hash_start, start, end, right_ext )
{
    if( ambig ) { done_ = true; return; }
    if( start >= end ) { done_ = true; return; }
    type_1_ = SErrType::D;
    alpha_ = alphabet.begin();
    hash_ = orig_hash_;
    prev_letter_ = 0xFF;

    if( right_ext_ ) {
        curr_pos_ = (hash_start_ + start_ > 0) ? start_ : 1;
        if( curr_pos_ >= end_ ) { done_ = true; return; }
        curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
        if( *alpha_ == prev_letter_ ) ++alpha_;
        InsertAndShiftRight();
    }
    else {
        curr_pos_ = (hash_start_ + end_ < q_->Len()) ? end_ : end_ - 1;
        if( curr_pos_ <= start_ ) { done_ = true; return; }
        curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ - 1 );
        if( *alpha_ == prev_letter_ ) ++alpha_;
        InsertAndShiftLeft();
    }
}

//------------------------------------------------------------------------------
bool CQueryData::CBlowUpIterator_D::Next( void )
{
    if( End() ) return false;
    hash_ = orig_hash_;
    ++alpha_;

    if( right_ext_ ) {
        if( alpha_ == alphabet.end() ) {
            if( ++curr_pos_ == end_ ) { done_ = true; return false; }
            alpha_ = alphabet.begin();
            prev_letter_ = curr_letter_;
            curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
            if( *alpha_ == prev_letter_ ) ++alpha_;
            InsertAndShiftRight();
        }
        else {
            if( *alpha_ == prev_letter_ ) {
                ++alpha_;
            }

            if( alpha_ == alphabet.end() ) {
                if( ++curr_pos_ == end_ ) { done_ = true; return false; }
                prev_letter_ = curr_letter_;
                curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ );
                alpha_ = alphabet.begin();
                if( *alpha_ == prev_letter_ ) ++alpha_;
            }

            InsertAndShiftRight();
        }
    }
    else {
        if( alpha_ == alphabet.end() ) {
            if( --curr_pos_ == start_ ) { done_ = true; return false; }
            prev_letter_ = curr_letter_;
            curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ - 1 );
            alpha_ = alphabet.begin();
            if( *alpha_ == prev_letter_ ) ++alpha_;
            InsertAndShiftLeft();
        }
        else {
            if( *alpha_ == prev_letter_ ) ++alpha_;

            if( alpha_ == alphabet.end() ) {
                if( --curr_pos_ == start_ ) { done_ = true; return false; }
                prev_letter_ = curr_letter_;
                curr_letter_ = GetLetter< CODING >( orig_hash_, curr_pos_ - 1 );
                alpha_ = alphabet.begin();
                if( *alpha_ == prev_letter_ ) ++alpha_;
            }

            InsertAndShiftLeft();
        }
    }

    return true;
}

//------------------------------------------------------------------------------
CQueryData::CBlowUpIterator::CBlowUpIterator( 
        const CQueryData & q, TSeqSize hash_start, 
        TSeqSize start, TSeqSize end, bool right_ext, bool ambig )
    : iter_m_( q, hash_start, start, end, right_ext, ambig ),
      iter_i_( q, hash_start, start, end, right_ext, ambig ),
      iter_d_( q, hash_start, start, end, right_ext, ambig ),
      curr_iter_( &iter_m_ )
{
    if( curr_iter_->End() ) curr_iter_ = &iter_i_;
    if( curr_iter_->End() ) curr_iter_ = &iter_d_;
}

//------------------------------------------------------------------------------
bool CQueryData::CBlowUpIterator::Next( void )
{
    if( End() ) return false;

    if( !curr_iter_->Next() ) {
        if( curr_iter_ == &iter_m_ ) {
            curr_iter_ = &iter_i_;

            if( curr_iter_->End() ) {
                curr_iter_ = &iter_d_;
                return !curr_iter_->End();
            }
            else return true;
        }
        else if( curr_iter_ == &iter_i_ ) {
            curr_iter_ = &iter_d_;
            return !curr_iter_->End();
        }
        else return false;
    }

    return true;
}

END_NS( srprism )
END_STD_SCOPES

