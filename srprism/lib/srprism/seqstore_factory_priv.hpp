/*  $Id: seqstore_factory_priv.hpp 351764 2012-02-01 14:07:34Z morgulis $
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
 * File Description: seqstore generation
 *
 */

#include "../common/trace.hpp"
#include "../seq/seqdef.hpp"
#include "seqstore_factory.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
template< typename data_t >
void CSeqStoreFactory::Append( const std::string & id, const data_t & seq_data )
{
    typedef SCodingTraits< SEQDATA_CODING > TSeqTraits;
    static const TLetter LETTER_MASK = SBitFieldTraits< 
        TLetter, TSeqTraits::LETTER_BITS >::MASK;
    typedef typename data_t::TCodingTraits TInputTraits;

    TLetter ambig_subst( 0 );
    bool in_ambig_run( false );
    M_TRACE( CTracer::INFO_LVL, 
             "adding sequence " << id << " of length " << seq_data.size );
    if( seq_data.size < min_seq_len_ ) min_seq_len_ = seq_data.size;
    TDBOrdId sid( id_map_.size() );
    id_map_.push_back( id );
    rev_id_map_.push_back( SRevIdMapEntry( id, sid ) );
    seq_info_.push_back( SSeqInfoEntry( sid, seq_data.size ) );

    {
        size_t dsz( (seq_data.size + WORD_LETTERS)/TSeqTraits::PACK_FACTOR );
        char * d( (char *)mem_mgr_.Allocate( dsz ) ),
             * dm( (char *)mem_mgr_.Allocate( dsz ) );
        std::fill( d, d + dsz, 0 );
        std::fill( dm, dm + dsz, 0 );
        seq_data_.push_back( (TWord *)d );
        mask_data_.push_back( (TWord *)dm );
    }

    for( TSeqSize k( 0 ); k < seq_data.size; k += WORD_LETTERS ) {
        TWord word( 0 ),
              mword( 0 );

        for( TSeqSize j( 0 ); j < WORD_LETTERS && j + k < seq_data.size; ++j ) {
            TLetter input_letter( seq_data.seq[k + j] );

            if( TInputTraits::NAMBIG[input_letter] ) {
                in_ambig_run = false;
                SetLetter< SEQDATA_CODING >( 
                        word, j, 
                        TInputTraits::template Recode< SEQDATA_CODING >(
                            input_letter ) );
            }
            else {
                if( !in_ambig_run ) {
                    ambig_map_.push_back(
                            SAmbigRunData( sid, k + j, ambig_data_.size() ) );
                    in_ambig_run = true;
                }
                
                ambig_data_.push_back( input_letter );
                SetLetter< SEQDATA_CODING >( word, j, ambig_subst );
                SetLetter< SEQDATA_CODING >( mword, j, 3 );
                ambig_subst = ((ambig_subst + 1)&LETTER_MASK);
                ++ambig_map_.rbegin()->len;
            }
        }

        (*seq_data_.rbegin())[k/WORD_LETTERS] = word;
        (*mask_data_.rbegin())[k/WORD_LETTERS] = mword;
    }
}

END_NS( srprism )
END_STD_SCOPES

