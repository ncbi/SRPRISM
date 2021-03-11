/*  $Id: seqstore.cpp 384196 2012-12-21 17:04:48Z morgulis $
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
 * File Description: part of the index containing packed subject sequence data
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <fstream>
#include <stdexcept>
#include <vector>

#include "../common/trace.hpp"
#include "../common/binfile.hpp"
#include "../common/util.hpp"
#include "seqstore.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
const char * CSeqStoreBase::DESC_SFX       = ".ssd";
const char * CSeqStoreBase::ID_MAP_SFX     = ".imp";
const char * CSeqStoreBase::POS_MAP_SFX    = ".pmp";
const char * CSeqStoreBase::SEQ_DATA_SFX   = ".ss";
const char * CSeqStoreBase::MASK_DATA_SFX  = ".mss";
const char * CSeqStoreBase::AMBIG_MAP_SFX  = ".amp";
const char * CSeqStoreBase::AMBIG_DATA_SFX = ".ssa";

//------------------------------------------------------------------------------
CSeqStore::CSeqStore( const std::string & basename, CMemoryManager & mem_mgr )
    : basename_( basename ), 
      n_seq_( 0 ), data_sz_( 0 ), ambig_map_sz_( 0 ), ambig_data_sz_( 0 ),
      mem_mgr_( mem_mgr ),
      ambig_map_( 0 ), ambig_data_( 0 ), seq_data_( 0 ),
      max_seq_overlap_( 0 )
{
}

//------------------------------------------------------------------------------
void CSeqStore::LoadHeader( void )
{
    CReadBinFile ins( basename_ + DESC_SFX );

    {
        Uint1 v( 0 );
        ins.Read( (char *)&v, 1, true );
        if( v != SS_VERSION ) M_THROW( CException, VALID, "version mismatch" );
        Uint1 junk[7];
        ins.Read( (char *)junk, 7, true );
    }

    Uint8 data( 0 );
    ins.Read( (char *)&data, 8, true );
    n_seq_ = (Uint4)data;
    ins.Read( (char *)&data, 8, true );
    data_sz_ = (Uint4)data;
    ins.Read( (char *)&data, 8, true );
    ambig_map_sz_ = (Uint4)data;
    ins.Read( (char *)&data, 8, true );
    ambig_data_sz_ = (Uint4)data;
    ins.Read( (char *)&segment_letters_, 8, true );
    ambig_mask_.resize( MaskUnits( segment_letters_ ), 0 );
    ins.Read( 
            (char *)&ambig_mask_[0], 
            sizeof( TMaskUnit )*MaskUnits( segment_letters_ ), 
            true );
    M_TRACE( CTracer::INFO_LVL, "sequence store header loaded" );
    M_TRACE( CTracer::INFO_LVL, "number of sequences: " << n_seq_ );
    M_TRACE( CTracer::INFO_LVL, 
             "sequence data size: " << data_sz_ << " letters" );
    M_TRACE( CTracer::INFO_LVL, 
             "ambiguity map size: " << ambig_map_sz_ << " regions" );
    M_TRACE( CTracer::INFO_LVL, 
             "ambiguity data size: " << ambig_data_sz_ << " letters" );
}

//------------------------------------------------------------------------------
namespace {
    struct SSegEnd
    {
        TDBOrdId sid;
        TSeqSize pos;
        TSeqSize osid;
        bool left;

        SSegEnd( 
                TDBOrdId s = 0, TDBOrdId o = 0, 
                TSeqSize p = 0, bool l = true )
            : sid( s ), pos( p ), osid( o ), left( l )
        {}

        friend bool operator<( const SSegEnd & l, const SSegEnd & r )
        { 
            if( l.sid == r.sid ) {
                if( l.pos == r.pos ) {
                    if( l.left && !r.left ) return false;
                    else if( r.left && !l.left ) return true;
                    else return true;
                }
                else return l.pos < r.pos;
            }
            return l.sid < r.sid;
        }
    };

    typedef std::vector< SSegEnd > TSegs;
}

void CSeqStore::ComputeSeqOverlap( void )
{
    TSegs segs;
    segs.reserve( 2*n_seq_ );

    for( size_t i( 0 ); i < n_seq_; ++i ) {
        SSeqMapEntry & s( seq_map_[i] );
        SSegEnd l( s.oid, i, s.ref_loc_start, true ),
                r( s.oid, i, s.ref_loc_end, false );
        segs.push_back( l );
        segs.push_back( r );
    }

    std::sort( segs.begin(), segs.end() );
    Uint4 overlap( 0 );

    for( TSegs::const_iterator i( segs.begin() ); i != segs.end(); ++i ) {
        if( i->left ) {
            ++overlap;
            if( max_seq_overlap_ < overlap ) max_seq_overlap_ = overlap;
        }
        else --overlap;
    }

    M_TRACE( CTracer::INFO_LVL, 
             "sequence overlap factor: " << max_seq_overlap_ );

    TALLists tmp_all;
    tmp_all.reserve( max_seq_overlap_ );

    for( TSegs::const_iterator i( segs.begin() ); i != segs.end(); ++i ) {
        if( i->left ) tmp_all.push_back( i->osid );
        else {
            TALLists::iterator f( std::find( 
                        tmp_all.begin(), tmp_all.end(), i->osid ) );
            if( f != tmp_all.end() ) tmp_all.erase( f );
        }

        TSegs::const_iterator j( i+1 );

        if( j == segs.end() || j->sid != i->sid || j->pos != i->pos ) {
            ref_segs_.push_back( 
                    SRefSegEnd( al_lists_.size(), i->sid, i->pos ) );

            for( TALLists::const_iterator j( tmp_all.begin() );
                    j != tmp_all.end(); ++j ) {
                al_lists_.push_back( *j );
            }
        }
    }

    M_TRACE( CTracer::INFO_LVL,
             "alternative loci overlap table size: " << al_lists_.size() );
    M_TRACE( CTracer::INFO_LVL,
             "alternative loci partition size: " << ref_segs_.size() );
}

//------------------------------------------------------------------------------
void CSeqStore::LoadSeqMap( void )
{
    seq_map_.resize( n_seq_ );
    CReadBinFile ins( basename_ + POS_MAP_SFX );

    for( size_t i( 0 ); i < n_seq_; ++i ) {
        SSeqMapEntry & s( seq_map_[i] );
        ins.Read( (char *)&s.oid, sizeof( TDBOrdId ), true );
        ins.Read( (char *)&s.seq_start, sizeof( TPos ), true );
        ins.Read( (char *)&s.body_start, sizeof( TPos ), true );
        ins.Read( (char *)&s.body_end, sizeof( TPos ), true );
        ins.Read( (char *)&s.seq_end, sizeof( TPos ), true );
        ins.Read( (char *)&s.ref_loc_start, sizeof( TSeqSize ), true );
        ins.Read( (char *)&s.ref_loc_end, sizeof( TSeqSize ), true );
    }

    M_TRACE( CTracer::INFO_LVL, "sequence map loaded" );
    ComputeSeqOverlap();
}

//------------------------------------------------------------------------------
void CSeqStore::LoadDynamicData( void )
{
    {
        ambig_map_ = 
            (SAmbigRun *)mem_mgr_.Allocate( ambig_map_sz_*sizeof( SAmbigRun ) );
        CReadBinFile ins( basename_ + AMBIG_MAP_SFX );
        ins.Read( (char *)ambig_map_, ambig_map_sz_*sizeof( SAmbigRun ), true );
        M_TRACE( CTracer::INFO_LVL, "ambiguity map loaded" );
    }

    size_t n_words( data_sz_/WORD_LETTERS );

    try { 
        seq_data_ = 
            (TWord *)mem_mgr_.Allocate( (n_words + 3)*sizeof( TWord ) );
        seq_data_[n_words + 1] = seq_data_[n_words + 2] = 0;
        *seq_data_++ = 0;
    }
    catch( ... ) {
        mem_mgr_.Free( (void *)ambig_map_ );
        ambig_map_ = 0;
        throw;
    }

    try {
        CReadBinFile ins( basename_ + SEQ_DATA_SFX );
        ins.Read( (char *)seq_data_, n_words*sizeof( TWord ), true );
        M_TRACE( CTracer::INFO_LVL, "sequence data loaded" );
    }
    catch( ... ) {
        mem_mgr_.Free( seq_data_ - 1 );
        seq_data_ = 0;
        mem_mgr_.Free( ambig_map_ );
        ambig_map_ = 0;
        throw;
    }

    try{ 
        rev_seq_data_ = 
            (TWord *)mem_mgr_.Allocate( (n_words + 3)*sizeof( TWord ) ); 
        rev_seq_data_[n_words + 1] = rev_seq_data_[n_words + 2] = 0;
        *rev_seq_data_++ = 0;
        
        for( size_t i( 0 ); i < n_words; ++i ) {
            ReverseComplement< SEQDATA_CODING >( 
                    rev_seq_data_[n_words - i - 1], seq_data_[i] );
        }

        M_TRACE( CTracer::INFO_LVL, "reverse sequence data generated" );
    }
    catch( ... ) {
        mem_mgr_.Free( seq_data_ - 1 );
        seq_data_ = 0;
        mem_mgr_.Free( ambig_map_ );
        ambig_map_ = 0;
        throw;
    }
}

//------------------------------------------------------------------------------
void CSeqStore::LoadAmbigData( void )
{
    if( ambig_data_ == 0 ) {
        ambig_data_ = (TLetter *)mem_mgr_.Allocate( 
                ambig_data_sz_*sizeof( TLetter ) );

        try {
            CReadBinFile ins( basename_ + AMBIG_DATA_SFX );
            ins.Read( 
                    (char *)ambig_data_, 
                    ambig_data_sz_*sizeof( TLetter ), 
                    true );
            M_TRACE( CTracer::INFO_LVL, "ambiguity data loaded" );
        }
        catch( ... ) {
            mem_mgr_.Free( ambig_data_ );
            ambig_data_ = 0;
            throw;
        }
    }
}

//------------------------------------------------------------------------------
void CSeqStore::UnloadAmbigData( void )
{
    if( ambig_data_ != 0 ) {
        mem_mgr_.Free( ambig_data_ );
        ambig_data_ = 0;
    }
}

//------------------------------------------------------------------------------
void CSeqStore::Load(void)
{
    if( seq_data_ == 0 ) {
        LoadHeader();
        LoadSeqMap();
        LoadDynamicData();
    }
}

//------------------------------------------------------------------------------
void CSeqStore::Unload( void )
{
    if( seq_data_ != 0 ) {
        mem_mgr_.Free( (void *)(rev_seq_data_ - 1));
        mem_mgr_.Free( (void *)(seq_data_ - 1));
        mem_mgr_.Free( (void *)ambig_map_ );
        rev_seq_data_ = seq_data_ = 0; ambig_map_ = 0;
    }
}

END_NS( srprism )
END_STD_SCOPES

