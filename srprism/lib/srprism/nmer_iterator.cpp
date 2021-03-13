/*  $Id: nmer_iterator.cpp 205515 2010-09-20 15:51:04Z morgulis $
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
 * File Description: sequence traversal by conseequtive n-mers
 *
 */

#include <ncbi_pch.hpp>

#include "../common/bits.hpp"
#include "nmer_iterator.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
bool CNMerIterator::NextLetter( void )
{
    static const size_t OFFSET_SHIFT = SBinLog< HASH_LEN >::VALUE;
    static const TSeqSize OFFSET_MASK = 
        SBitFieldTraits< TSeqSize, OFFSET_SHIFT >::MASK;
    if( End() || left_ < HASH_LEN ) { end_ = true; return false; }
    curr_prefix_ = GetWord< SEQDATA_CODING >( fw_ptr_, fw_off_ );
    ReverseComplement< SEQDATA_CODING >( reverse_prefix_, curr_prefix_ );
    len_ = (ss_.NAmbigs( curr_pos_, curr_pos_ + HASH_LEN ) > 0 ) ? 0 : HASH_LEN;
    // len_ = HASH_LEN;
    // len_ = (ss_.NAmbigs( curr_pos_, curr_pos_ + HASH_LEN ) > 15 ) ? 0 : HASH_LEN;
    ++curr_pos_; --left_; ++fw_off_; fw_ptr_ += (fw_off_>>OFFSET_SHIFT);
    fw_off_ &= OFFSET_MASK;
    return true;
}

//------------------------------------------------------------------------------
CNMerIterator::CNMerIterator( const CSeqStore & ss, TDBOrdId seq_n )
    : ss_( ss ), len_( 0 ), left_( 0 ), 
      fw_off_( 0 ), curr_pos_( 0 ), fw_ptr_( 0 ),
      curr_prefix_( 0 ), reverse_prefix_( 0 ), end_( false )
{ 
    curr_pos_ = ss_.EncodePos( std::make_pair( seq_n, 0 ) );
    left_ = ss_.FwTailLen( curr_pos_ );
    std::pair< const TWord *, TSeqSize > t( ss_.FwDataPtr( curr_pos_ ) );
    fw_ptr_ = t.first; fw_off_ = t.second;
    end_ = (left_ < HASH_LEN);
    Next(); 
}


END_NS( srprism )
END_STD_SCOPES

