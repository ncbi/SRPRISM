/*  $Id: seqstore_factory.cpp 431273 2014-04-02 17:10:44Z morgulis $
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

#include <ncbi_pch.hpp>

#include <algorithm>

#include "../common/bits.hpp"
#include "../common/util.hpp"
#include "../common/textfile.hpp"
#include "seqstore_factory.hpp"

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
CSeqStoreFactory::CSeqStoreFactory( 
        size_t max_mem, 
        const std::string & base_name, 
        const std::string & alt_loc_spec_name,
        size_t segment_letters,
        size_t al_extend )
    : mem_mgr_( max_mem ), base_name_( base_name ), 
      alt_loc_spec_name_( alt_loc_spec_name ),
      ss_outs_( 0 ), curr_pos_( 0 ),
      ambig_map_size_( 0 ), ambig_data_size_( 0 ),
      segment_letters_( segment_letters ),
      min_seq_len_( common::SIntTraits< size_t >::MAX ),
      al_extend_( (TSeqSize)al_extend )
{ 
    ss_outs_.reset( new CWriteBinFile( base_name + SEQ_DATA_SFX) );
    seqmap_outs_.reset( new CWriteBinFile( base_name + POS_MAP_SFX ) );
    ambig_map_outs_.reset( new CWriteBinFile( base_name + AMBIG_MAP_SFX ) );
    ambig_data_outs_.reset( new CWriteBinFile( base_name + AMBIG_DATA_SFX ) );
    idmap_outs_ = CWriteTextFile::MakeWriteTextFile( base_name + ID_MAP_SFX );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::ComputeSeqLetters( void )
{
    static const size_t MIN_SEG_LETTERS = 32;
    static const size_t MAX_SEG_LETTERS = 8192;

    if( segment_letters_ == 0 ) {
        segment_letters_ = MIN_SEG_LETTERS;

        while( segment_letters_ <= MAX_SEG_LETTERS && 
                16*segment_letters_ < min_seq_len_  ) {
            segment_letters_ *= 2;
        }
    }

    M_TRACE( CTracer::INFO_LVL, 
             "using segment size of " << segment_letters_ << " bases" );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::SaveHeader( void )
{
    M_TRACE( CTracer::INFO_LVL, "saving header" );
    CWriteBinFile hdr_outs( base_name_ + DESC_SFX );

    {
        Uint1 ver( (Uint1)SS_VERSION );
        hdr_outs.Write( (const char *)&ver, 1 );
    }

    {
        char a[7] = {0, 0, 0, 0, 0, 0, 0};
        hdr_outs.Write( (const char *)a, 7 );
    }

    Uint8 data( (Uint8)id_map_.size() );
    hdr_outs.Write( (const char *)&data, 8 );
    data = (Uint8)curr_pos_;
    hdr_outs.Write( (const char *)&data, 8 );
    data = (Uint8)ambig_map_size_;
    hdr_outs.Write( (const char *)&data, 8 );
    data = (Uint8)ambig_data_size_;
    hdr_outs.Write( (const char *)&data, 8 );
    hdr_outs.Write( (const char *)&segment_letters_, 8 );
    hdr_outs.Write( 
            (const char *)&ambig_mask_[0], 
            MaskUnits( segment_letters_ )*sizeof( TMaskUnit ) );
}

//------------------------------------------------------------------------------
std::string CSeqStoreFactory::GetWord( 
        const std::string & line, std::string::size_type & pos )
{
    std::string res;
    pos = line.find_first_not_of( " \t", pos );

    if( pos != std::string::npos ) {
        std::string::size_type pos1( line.find_first_of( " \t", pos ) );
        res = line.substr( pos, pos1 - pos );
        pos = pos1 + 1;
    }

    return res;
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::ParseAltLocLine( const std::string & line )
{
    std::string::size_type pos( 0 ); GetWord( line, pos );
    std::string word( GetWord( line, pos ) );
    TDBOrdId oid, ref_oid;

    // check for presence of subject id in the input
    //
    {
        TRevIdMap::const_iterator i( 
                std::lower_bound( rev_id_map_.begin(), rev_id_map_.end(), 
                SRevIdMapEntry( word, 0 ) ) );

        if( i == rev_id_map_.end() || i->id != word ) {
            M_THROW( CException, ALTLOC, 
                     "alternative locus id " << word << 
                     " is not present in the input" );
        }

        oid = i->oid;
    }

    word = GetWord( line, pos );

    // check for presence of reference subject id in the input
    //
    {
        TRevIdMap::const_iterator i( 
                std::lower_bound( rev_id_map_.begin(), rev_id_map_.end(), 
                SRevIdMapEntry( word, 0 ) ) );

        if( i == rev_id_map_.end() || i->id != word ) {
            M_THROW( CException, ALTLOC, 
                     "reference id " << word << 
                     " is not present in the input" );
        }

        ref_oid = i->oid;
    }

    word = GetWord( line, pos );
    bool flip( false );

    if( word == "-" ) flip = true;
    else if( word != "+" ) {
        M_THROW( CException, ALTLOC,
                 "wrong value of strand field (must be '+' or '-') in " <<
                 line );
    }

    word = GetWord( line, pos );
    TSeqSize ref_loc_start( atoi( word.c_str() ) - 1 );
    word = GetWord( line, pos );
    bool left_fuzzy( word == "y" );

    if( word != "y" && word != "n" ) {
        M_THROW( CException, ALTLOC,
                 "wrong value of left fuzzy indicator field "
                 "(must be 'y' or 'n') in " <<
                 line );
    }

    word = GetWord( line, pos );
    TSeqSize ref_loc_end( atoi( word.c_str() ) );
    word = GetWord( line, pos );
    bool right_fuzzy( word == "y" );

    if( word != "y" && word != "n" ) {
        M_THROW( CException, ALTLOC,
                 "wrong value of right fuzzy indicator field "
                 "(must be 'y' or 'n') in " <<
                 line );
    }

    word = GetWord( line, pos );
    TSeqSize loc_start( atoi( word.c_str() ) - 1 );
    word = GetWord( line, pos );
    TSeqSize loc_end( atoi( word.c_str() ) );

    if( loc_start >= loc_end ) {
        M_THROW( CException, ALTLOC,
                 "start coordinate of the alternate locus " << loc_start << 
                 " must be less or equal to the end coordinate " << loc_end << 
                 " in line " << line );
    }

    if( ref_loc_start > ref_loc_end ) {
        M_THROW( CException, ALTLOC,
                 "start coordinate of the reference " << ref_loc_start << 
                 " must be less or equal to the end coordinate of the "
                 "reference " << ref_loc_end << " in line " << line );
    }

    seq_info_[oid] = SSeqInfoEntry( 
            oid, ref_oid, loc_start, loc_end, 
            seq_info_[oid].alt_loc_real_end, ref_loc_start, ref_loc_end,
            left_fuzzy, right_fuzzy, flip );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::SetUpSeqInfo( void )
{
    if( alt_loc_spec_name_.empty() ) return;
    M_TRACE( CTracer::INFO_LVL, 
             "reading alternate loci specification from " << 
             alt_loc_spec_name_ );
    std::auto_ptr< CReadTextFile > al_is( 
            CReadTextFile::MakeReadTextFile( alt_loc_spec_name_ ) );

    while( !al_is->Eof() ) {
        std::string al_line( al_is->GetLine() );
        if( al_line.empty() || al_line[0] == '#' ) continue;
        ParseAltLocLine( al_line );
    }

    std::sort( seq_info_.begin(), seq_info_.end() );
    M_TRACE( CTracer::INFO_LVL, "alternate loci descriptor table created" );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::StoreSeqSeg( 
        const TWord * seg_data, TSeqSize seg_start, TSeqSize seg_len, 
        TWord * dest, TSeqSize dest_start )
{
    if( seg_len > 0 ) {
        for( TSeqSize j( 0 ); j < seg_len; ++j ) {
            TLetter l( GetStreamLetter< SEQDATA_CODING >( 
                        seg_data, seg_start + j ) );
            SetStreamLetter< SEQDATA_CODING >( dest, dest_start + j, l );
        }
    }
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::StoreAmbigInfoReverse( 
        TDBOrdId oid, TSeqSize start_pos, TSeqSize len, TSeqSize offset,
        std::vector< SAmbigRun > & amap, std::vector< TLetter > & adata )
{
    typedef SCodingTraits< OUTPUT_CODING > TOutCodeTraits;
    if( len == 0 || ambig_map_.empty() ) return;
    SAmbigRunData p( oid, start_pos + len, 0 );
    TAmbigMap::const_iterator i(
            std::upper_bound( ambig_map_.begin(), ambig_map_.end(), p ) );
    if( i == ambig_map_.begin() ) return;

    do {
        --i;

        if( i->sid < oid || i->pos + i->len < start_pos ) break;
        TPos s( std::max( start_pos, i->pos ) ),
             e( std::min( start_pos + len, i->pos + i->len ) );
        SRPRISM_ASSERT( s < e );
        SAmbigRun t( curr_pos_ + offset + (start_pos + len - e),
                     ambig_data_size_ + adata.size() );
        t.len = (e - s);
        amap.push_back( t );
        TSeqSize off( start_pos - i->pos );

        for( TSeqSize j( t.len ); j > 0; --j ) {
            TLetter l( ambig_data_[i->offset + off + (j-1)] );
            adata.push_back( TOutCodeTraits::RC[l] );
        }
    }
    while( i != ambig_map_.begin() );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::StoreAmbigInfo( 
        TDBOrdId oid, TSeqSize start_pos, TSeqSize len, TSeqSize offset,
        std::vector< SAmbigRun > & amap, std::vector< TLetter > & adata )
{
    if( len > 0 ) {
        SAmbigRunData p( oid, start_pos, 0 );
        TAmbigMap::const_iterator i( 
                std::lower_bound( ambig_map_.begin(), ambig_map_.end(), p ) );
        
        if( i != ambig_map_.begin() ) {
            --i;
            
            if( i->sid == p.sid && i->pos + i->len > p.pos ) {
                SAmbigRun t( 
                        curr_pos_ + offset, ambig_data_size_ + adata.size() );
                t.len = i->pos + i->len - start_pos;
                amap.push_back( t );

                for( size_t j( 0 ); j < t.len; ++j ) {
                    adata.push_back( 
                            ambig_data_[i->offset + (start_pos - i->pos) + j] );
                }
            }

            ++i;
        }

        for( ; i != ambig_map_.end() && 
               i->sid == oid && 
               i->pos < start_pos + len; ++i ) {
            TSeqSize end( std::min( i->pos + i->len, start_pos + len ) );
            SAmbigRun t( 
                    curr_pos_ + offset + (i->pos - start_pos), 
                    ambig_data_size_ + adata.size() );
            t.len = end - i->pos;
            amap.push_back( t );

            for( size_t j( 0 ); j < t.len; ++j ) {
                adata.push_back( ambig_data_[i->offset + j] );
            }
        }
    }
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::SaveSeqData( TDBOrdId oid )
{
    SSeqInfoEntry & si( seq_info_[oid] ), & ref_si( seq_info_[si.ref_oid] );
    bool left_fuzzy( si.fuzzy_left ), right_fuzzy( si.fuzzy_right );

    if( si.alt_loc_start != 0 ) {
        M_TRACE( CTracer::WARNING_LVL,
                 "start of alternate locus " << oid << 
                 " is moved to the start of the actual sequence" );
        si.alt_loc_start = 0;
    }

    if( si.alt_loc_end != si.alt_loc_real_end ) {
        M_TRACE( CTracer::WARNING_LVL,
                 "end of alternate locus " << oid << 
                 " is moved to the end of the actual sequence" );
        si.alt_loc_end = si.alt_loc_real_end;
    }

    TSeqSize pfx_len( std::min( al_extend_, si.ref_loc_start ) ),
             sfx_len( std::min( al_extend_,
                                ref_si.alt_loc_end 
                                    - ref_si.alt_loc_start 
                                    - si.ref_loc_end ) );

    if( si.flip ) { 
        std::swap( pfx_len, sfx_len );
        std::swap( left_fuzzy, right_fuzzy );
    }

    if( left_fuzzy ) pfx_len = 0;
    if( right_fuzzy ) sfx_len = 0;

    TSeqSize pfx_start_pos, sfx_start_pos;

    if( si.flip ) {
        sfx_start_pos = si.ref_loc_start - sfx_len;
        pfx_start_pos = si.ref_loc_end;
    }
    else {
        pfx_start_pos = si.ref_loc_start - pfx_len;
        sfx_start_pos = si.ref_loc_end;
    }

    TSeqSize body_len( si.alt_loc_end - si.alt_loc_start );
    TSeqSize full_len( pfx_len + body_len + sfx_len );
    size_t n_segs( 1 + (full_len + WORD_LETTERS - 1)/segment_letters_ );
    size_t dsz( 
            1 + (n_segs*segment_letters_)/SCodingTraits< 
                SEQDATA_CODING >::PACK_FACTOR );
    TWord * d( (TWord *)mem_mgr_.Allocate( dsz ) );
    std::fill( (char *)d, (char*)d + dsz, 0 );

    // generate extended sequence data
    //
    TDBOrdId pfx_id( left_fuzzy ? si.oid : si.ref_oid ),
             sfx_id( right_fuzzy ? si.oid : si.ref_oid );

    if( si.flip ) {
        TWord * t( new TWord[pfx_len/WORD_LETTERS + 1] );
        ReverseComplement< SEQDATA_CODING >( 
                t, seq_data_[pfx_id], pfx_start_pos, pfx_len );
        StoreSeqSeg( t, 0, pfx_len, d, 0 );
        delete[] t;
    }
    else StoreSeqSeg( seq_data_[pfx_id], pfx_start_pos, pfx_len, d, 0 );

    StoreSeqSeg( seq_data_[si.oid], si.alt_loc_start, body_len, d, pfx_len );

    if( si.flip ) {
        TWord * t( new TWord[sfx_len/WORD_LETTERS + 1] );
        ReverseComplement< SEQDATA_CODING >(
                t, seq_data_[sfx_id], sfx_start_pos, sfx_len );
        StoreSeqSeg( t, 0, sfx_len, d, pfx_len + body_len );
        delete[] t;
    }
    else {
        StoreSeqSeg( 
                seq_data_[sfx_id], sfx_start_pos, sfx_len, 
                d, pfx_len + body_len );
    }

    // create ambiguity information
    //
    std::vector< SAmbigRun > amap;
    std::vector< TLetter > adata;

    if( si.flip ) {
        StoreAmbigInfoReverse( pfx_id, pfx_start_pos, pfx_len, 0, amap, adata );
        StoreAmbigInfo( si.oid, si.alt_loc_start, 
                        body_len, pfx_len, amap, adata );
        StoreAmbigInfoReverse( sfx_id, sfx_start_pos, 
                               sfx_len, pfx_len + body_len, amap, adata );
    }
    else {
        StoreAmbigInfo( pfx_id, pfx_start_pos, pfx_len, 0, amap, adata );
        StoreAmbigInfo( si.oid, si.alt_loc_start, 
                        body_len, pfx_len, amap, adata );
        StoreAmbigInfo( sfx_id, si.ref_loc_end, 
                        sfx_len, pfx_len + body_len, amap, adata );
    }

    // update ambiguous segment mask
    //
    if( !amap.empty() ) {
        for( size_t i( 0 ); i < n_segs; ++i ) {
            SAmbigRun t( curr_pos_ + i*segment_letters_, 0 );
            std::vector< SAmbigRun >::const_iterator j(
                    std::lower_bound( amap.begin(), amap.end(), t ) );

            if( (j != amap.end() && j->pos < t.pos + segment_letters_ ) ||
                (j != amap.begin() && (j-1)->pos + (j-1)->len >= t.pos) ) {
                size_t seg_idx( curr_pos_/segment_letters_ + i );
                size_t unit( seg_idx/MASK_UNIT_BITS ),
                       bit_idx( seg_idx%MASK_UNIT_BITS );
                AssignBit( bit_idx, ambig_mask_[unit], true );
            }
        }
    }

    // append the computed data to the relevant files
    //
    M_TRACE( CTracer::INFO_LVL, 
             "sequence id: " << id_map_[si.oid] <<
             "; reference id: " << id_map_[si.ref_oid] );
    idmap_outs_->LineOut( id_map_[si.oid] );

    M_TRACE( CTracer::INFO_LVL,
             "ambiguous regions: " << amap.size() <<
             "; ambiguous bases: " << adata.size() );

    {
        Uint4 sz( (Uint4)amap.size() );

        if( sz > 0 ) {
            ambig_map_outs_->Write( 
                    (const char *)&amap[0], sz*sizeof( SAmbigRun ) );
        }
    }

    {
        Uint4 sz( (Uint4)adata.size() );

        if( sz > 0 ) {
            ambig_data_outs_->Write( 
                    (const char *)&adata[0], sz*sizeof( TLetter ) );
        }
    }

    {
        TPos p1( curr_pos_ ), p2( p1 + pfx_len ), 
             p3( p2 + body_len ), p4( p3 + sfx_len );
        M_TRACE( CTracer::INFO_LVL,
                 "sequence start: " << p1 << "; body start: " << p2 <<
                 "; body end: " << p3 << "; sequence end: " << p4 );
        seqmap_outs_->Write( (const char *)&si.ref_oid, sizeof( TDBOrdId ) );
        seqmap_outs_->Write( (const char *)&p1, sizeof( TPos ) );
        seqmap_outs_->Write( (const char *)&p2, sizeof( TPos ) );
        seqmap_outs_->Write( (const char *)&p3, sizeof( TPos ) );
        seqmap_outs_->Write( (const char *)&p4, sizeof( TPos ) );
        seqmap_outs_->Write( 
                (const char *)&si.ref_loc_start, sizeof( TSeqSize ) );
        seqmap_outs_->Write( 
                (const char *)&si.ref_loc_end, sizeof( TSeqSize ) );
    }

    ss_outs_->Write( 
            (const char *)d, 
            n_segs*SegmentWords( segment_letters_ )*sizeof( TWord ) );

    if( DATA_SIZE_LIMIT <= (size_t)curr_pos_ + n_segs*segment_letters_ ) {
        M_THROW( CException, MEMORY, 
                    "seqstore overflow: please, consider creating "
                    "several smaller databases" );
    }

    curr_pos_ += n_segs*segment_letters_;
    ambig_map_size_ += amap.size();
    ambig_data_size_ += adata.size();
    mem_mgr_.Free( d );
}

//------------------------------------------------------------------------------
void CSeqStoreFactory::Save( void )
{
    ComputeSeqLetters();
    ambig_mask_.resize( MaskUnits( segment_letters_ ), 0 );
    idmap_outs_->LineOut( id_map_.size() );
    std::sort( rev_id_map_.begin(), rev_id_map_.end() );
    M_TRACE( CTracer::INFO_LVL, "reverse id map generated" );
    SetUpSeqInfo();
    for( TDBOrdId i( 0 ); i < id_map_.size(); ++i ) SaveSeqData( i );
    SaveHeader();
}

END_NS( srprism )
END_STD_SCOPES

