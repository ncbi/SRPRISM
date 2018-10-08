/*  $Id: mkidx_pass.cpp 426095 2014-02-05 20:58:39Z morgulis $
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
 * File Description: one pass of index creation
 *
 */

#include <ncbi_pch.hpp>

#include <cassert>
#include <algorithm>

#include "../common/trace.hpp"
#include "../common/bits.hpp"
#include "../seq/seqdef.hpp"
#include "nmer_iterator.hpp"
#include "mkidx_pass.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
namespace {
    size_t NBytes( size_t v )
    {
        if( v == 0 ) return 0;
        else if( v <= 0xFF ) return 1;
        else if( v <= 0xFFFF ) return 2;
        else if( v <= 0xFFFFFFFF ) return 4;
        else return 8;
    }

    bool CheckPalindrome( TPrefix prefix )
    {
        TPrefix rc_prefix( 0 );
        ReverseComplement< SEQDATA_CODING >( rc_prefix, prefix );
        return (prefix == rc_prefix);
    }

    template< typename int_t >
    inline void VWrite( CWriteBinFile & os, size_t val )
    {
        // TODO: add range check.
        int_t v( (int_t)val );
        os.Write( (const char *)&v, sizeof( v ) );
    }

    inline void VWrite( CWriteBinFile & os, size_t bytes, size_t val )
    {
        switch( bytes ) {
            case 0: break;
            case 1: VWrite< Uint1 >( os, val ); break;
            case 2: VWrite< Uint2 >( os, val ); break;
            case 4: VWrite< Uint4 >( os, val ); break;
            case 8: VWrite< Uint8 >( os, val ); break;
            default: SRPRISM_ASSERT( false );
        }
    }
}

//------------------------------------------------------------------------------
std::pair< size_t, size_t > CMkIdxPass::SplitByStrand( 
        SPosEntry * entries, size_t start, size_t end )
{
    size_t f( 0 ), r( end - start );
    for( ; start != end && entries[start].strand == STRAND_FW; ++start ) ++f;
    r -= f;
    return std::make_pair( f, r );
}

//------------------------------------------------------------------------------
CMkIdxPass::CMkIdxPass( 
        const size_t * counts_table, const CSeqStore & seq_store, 
        void * free_space, size_t free_space_size, 
        size_t & start_idx, CWriteBinFile & map_file, CWriteBinFile & rmap_file,
        CWriteBinFile & idx_file )
    : seq_store_( seq_store ), map_file_( map_file ), rmap_file_( rmap_file ),
      idx_file_( idx_file ), total_entries_( 0 ), rmap_size_( 0 )
{
    M_TRACE( CTracer::INFO_LVL, 
             "start prefix: " << std::hex << start_idx << std::dec );
    size_t total( 0 );
    size_t start_idx_orig( start_idx );

    while( start_idx < NUM_HASH_KEYS ) {
        size_t hk_sz( counts_table[start_idx]*sizeof( SPosEntry ) );
        if( total + hk_sz >= free_space_size ) break;
        total += hk_sz;
        total_entries_ += counts_table[start_idx];
        ++start_idx;
    }

    M_TRACE( CTracer::INFO_LVL,
             "end prefix: " << std::hex << start_idx << std::dec );
    size_t start_prefix( start_idx_orig<<SHIFT ),
           end_prefix( start_idx<<SHIFT );
    SPosEntry * pos_data( (SPosEntry *)free_space ), 
              * pos_data_end( pos_data );
    M_TRACE( CTracer::INFO_LVL, "collecting position data" );

    for( size_t seq_idx = 0; seq_idx < seq_store_.NSeq(); ++seq_idx ) {
        for( CNMerIterator nmer_iter( seq_store, seq_idx );
                !nmer_iter.End(); nmer_iter.Next() ) {
            TPrefix prefix( nmer_iter.Prefix() );

            if( prefix >= start_prefix && prefix < end_prefix ) {
                SRPRISM_ASSERT( (size_t)(pos_data_end - pos_data) < total_entries_ );
                SPosEntry & e( *pos_data_end );
                e.prefix = prefix; 
                e.strand = nmer_iter.Strand();
                e.pos = nmer_iter.Pos() - HASH_LEN;
                ++pos_data_end;
            }
        }
    }

    SRPRISM_ASSERT( total_entries_ == (size_t)(pos_data_end - pos_data) );
    M_TRACE( CTracer::INFO_LVL, "loaded " << total_entries_ << " positions" );
    std::sort( pos_data, pos_data_end );
    entries_ = pos_data;
    start_idx_ = start_idx_orig;
    end_idx_ = start_idx;
}

//------------------------------------------------------------------------------
size_t CMkIdxPass::EstimateNormal( size_t start, size_t end )
{
    std::pair< size_t, size_t > strand_counts( 
            SplitByStrand( entries_, start, end ) );
    size_t & f( strand_counts.first ), & r( strand_counts.second );
    return NBytes( f ) + NBytes( r ) + (f + r)*4;
}

//------------------------------------------------------------------------------
size_t CMkIdxPass::CountExtensions(
        SPosEntry * entries, size_t start, size_t end ) const
{
    size_t n_ext( 0 );

    for( TExtensionIterator i( entries, start, end ); !i.End(); i.Next() ) {
        ++n_ext;
    }
    
    return n_ext;
}

//------------------------------------------------------------------------------
size_t CMkIdxPass::EstimateExtensions( 
        SPosEntry * entries, size_t start, size_t end ) const
{
    size_t n_ext( CountExtensions( entries, start, end ) );
    return 8 + n_ext*16 + (end - start)*4;
}

//------------------------------------------------------------------------------
void CMkIdxPass::SetUpExtensions( 
        SPosEntry * entries, size_t start, size_t end, TStrand s )
{
    for( size_t i = start; i < end; ++i ) {
        SPosEntry & entry( entries[i] );
        TStrand es( CombineStrands( s, entry.strand ) );
        std::pair< const TWord *, TSeqSize > t;

        if( es == STRAND_FW ) t = seq_store_.FwDataPtr( entry.pos + HASH_LEN );
        else t = seq_store_.RvDataPtr( entry.pos );

        entry.extension = GetWord< SEQDATA_CODING >( t.first, t.second );
    }

    std::sort( entries + start, entries + end );
}

//------------------------------------------------------------------------------
void CMkIdxPass::SetUpPalindromeData( 
        std::vector< SPosEntry > & entries, size_t start, size_t end )
{
    std::copy( entries_ + start, entries_ + end, entries.begin() );
    std::copy( 
            entries_ + start, entries_ + end, 
            entries.begin() + (end - start) );
        
    for( size_t i = end - start; i < entries.size(); ++i ) {
        entries[i].strand = STRAND_RV;
    }

    SetUpExtensions( &entries[0], 0, entries.size(), STRAND_FW );
}

//------------------------------------------------------------------------------
size_t CMkIdxPass::EstimateSpecial( size_t start, size_t end )
{
    size_t res( 0 );

    if( CheckPalindrome( entries_[start].prefix ) ) {
        std::vector< SPosEntry > entries( 2*(end - start) );
        SetUpPalindromeData( entries, start, end );
        res += 8 + EstimateExtensions( &entries[0], 0, entries.size() );
    }
    else {
        TPrefix prefix( entries_[start].prefix );
        SetUpExtensions( entries_, start, end, STRAND_FW );
        res += EstimateExtensions( entries_, start, end );
        SetUpExtensions( entries_, start, end, STRAND_RV );
        res += EstimateExtensions( entries_, start, end );
        for( size_t i = start; i < end; ++i ) entries_[i].prefix = prefix;
    }

    return res;
}

//------------------------------------------------------------------------------
void CMkIdxPass::UpdateRMap( size_t start, size_t end )
{
    SRMapEntry rmap_entry;
    rmap_entry.prefix = (TPrefix)entries_[start].prefix;
    std::pair< size_t, size_t > strand_counts( 
            SplitByStrand( entries_, start, end ) );
    rmap_entry.flog = (Uint1)BinLog( 1 + strand_counts.first );
    rmap_entry.rlog = (Uint1)BinLog( 1 + strand_counts.second );

    if( std::max( rmap_entry.flog, rmap_entry.rlog ) >= NPOS_LOG_START ) {
        rmap_file_.Write( 
                (const char *)&rmap_entry.prefix, sizeof( TPrefix ) );
        rmap_file_.Write(
                (const char *)&rmap_entry.flog, sizeof( Uint2 ) );
        rmap_file_.Write(
                (const char *)&rmap_entry.rlog, sizeof( Uint2 ) );
        ++rmap_size_;
    }
}

//------------------------------------------------------------------------------
size_t CMkIdxPass::Estimate( size_t start_entry, size_t end_entry )
{
    size_t res( 0 );

    for( CPrefixIterator pit( entries_, start_entry, end_entry );
            !pit.End(); pit.Next() ) {
        res += 2; // descriptor length
        size_t se( pit.StartEntry() ), ee( pit.EndEntry() );
        UpdateRMap( se, ee );
        
        if( ee - se >= FORMAT_THRESHOLD ) res += EstimateSpecial( se, ee );
        else res += EstimateNormal( se, ee );
    }

    return res;
}

//------------------------------------------------------------------------------
void CMkIdxPass::FlushSpecialForStrand( 
        SPosEntry * entries, size_t start, size_t end )
{
    size_t n_ext( CountExtensions( entries, start, end ) );
    SRPRISM_ASSERT( n_ext < 0x100000000ULL );

    {
        Uint4 t( (Uint4)n_ext );
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        DBG_TRACE( "IDXOUT: next: " << t );
        SRPRISM_ASSERT( end - start < 0x100000000ULL );
        t = (Uint4)(end - start);
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        DBG_TRACE( "IDXOUT: npos: " << t );
    }

    for( TExtensionIterator ei( entries, start, end ); !ei.End(); ei.Next() ) {
        size_t se( ei.StartEntry() ), ee( ei.EndEntry() );
        std::pair< size_t, size_t > s_counts( 
                SplitByStrand( entries, se, ee ) );
        SPosEntry & e( entries[se] );
        idx_file_.Write( (const char *)&e.extension, sizeof( e.extension ) );
        DBG_TRACE( "IDXOUT: ext: " << e.extension );
        SRPRISM_ASSERT( s_counts.first  < 0x100000000ULL );
        SRPRISM_ASSERT( s_counts.second < 0x100000000ULL );
        Uint4 t( (Uint4)s_counts.first );
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        DBG_TRACE( "IDXOUT: fnpos: " << t );
        t = (Uint4)s_counts.second;
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        DBG_TRACE( "IDXOUT: rnpos: " << t );
        SRPRISM_ASSERT( se < 0x100000000ULL );
        t = (Uint4)(se - start);
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        DBG_TRACE( "IDXOUT: idx: " << t );
    }

    for( ; start != end; ++start ) {
        idx_file_.Write( 
                (const char *)&entries[start].pos,
                sizeof( entries[start].pos ) );
        DBG_TRACE( "IDXOUT: pos: " << entries[start].pos );
    }
}

//------------------------------------------------------------------------------
void CMkIdxPass::FlushSpecial( size_t start, size_t end )
{
    DBG_TRACE( 
            "IDXOUT: prefix: " << std::hex << 
            entries_[start].prefix << std::dec );

    if( CheckPalindrome( entries_[start].prefix ) ) {
        DBG_TRACE( "IDXOUT: palindrome" );
        std::vector< SPosEntry > entries( 2*(end - start) );
        SetUpPalindromeData( entries, start, end );
        FlushSpecialForStrand( &entries[0], 0, entries.size() );
        Uint4 t( 0 );
        idx_file_.Write( (const char *)&t, sizeof( t ) );
        idx_file_.Write( (const char *)&t, sizeof( t ) );
    }
    else {
        SetUpExtensions( entries_, start, end, STRAND_FW );
        DBG_TRACE( "IDXOUT: strand 0" );
        FlushSpecialForStrand( entries_, start, end );
        SetUpExtensions( entries_, start, end, STRAND_RV );
        DBG_TRACE( "IDXOUT: strand 1" );
        FlushSpecialForStrand( entries_, start, end );
    }
}

//------------------------------------------------------------------------------
void CMkIdxPass::FlushNormal( size_t start, size_t end, Uint2 descr )
{
    DBG_TRACE( 
            "IDXOUT: prefix: " << std::hex << 
            entries_[start].prefix << std::dec );

    std::pair< size_t, size_t > strand_counts( 
            SplitByStrand( entries_, start, end ) );
    size_t & f( strand_counts.first ), & r( strand_counts.second );
    size_t f_bytes( NBytes( f ) );
    size_t r_bytes( NBytes( r ) );
    SetField< R_BYTES_START_BIT, R_BYTES_END_BIT >( descr, (Uint2)r_bytes );
    SetField< F_BYTES_START_BIT, F_BYTES_END_BIT >( descr, (Uint2)f_bytes );
    idx_file_.Write( (const char *)&descr, sizeof( descr ) );
    DBG_TRACE( "IDXOUT: descriptor: " << std::hex << descr << std::dec );
    VWrite( idx_file_, f_bytes, f );
    DBG_TRACE( "IDXOUT: fnpos: " << f );
    VWrite( idx_file_, r_bytes, r );
    DBG_TRACE( "IDXOUT: rnpos: " << r );

    for( ; start != end; ++start ) {
        idx_file_.Write( 
                (const char *)&entries_[start].pos,
                sizeof( entries_[start].pos ) );
        DBG_TRACE( "IDXOUT: pos: " << entries_[start].pos );
    }
}

//------------------------------------------------------------------------------
void CMkIdxPass::FlushMapPrefix( size_t start_entry, size_t end_entry )
{
    static const size_t SFX_MASK = SBitFieldTraits< 
            size_t, LETTER_BITS*(PREFIX_LEN - MAP_PREFIX_LEN) >::MASK;

    for( CPrefixIterator pit( entries_, start_entry, end_entry );
            !pit.End(); pit.Next() ) {
        size_t se( pit.StartEntry() ), ee( pit.EndEntry() );
        Uint2 descr( ((entries_[se].prefix)&SFX_MASK)<<SFX_START_BIT );

        if( ee - se  >= FORMAT_THRESHOLD ) {
            idx_file_.Write( (const char *)&descr, sizeof( descr ) );
            DBG_TRACE( 
                    "IDXOUT: descriptor: " << std::hex << descr << std::dec );
            FlushSpecial( se, ee );
        }
        else FlushNormal( se, ee, descr );
    }
}

//------------------------------------------------------------------------------
void CMkIdxPass::Run( void )
{
    CMapPrefixIterator map_prefixes( entries_, total_entries_ );
    size_t curr_prefix( start_idx_<<(SHIFT - CMapPrefixIterator::MASK_BITS) );
    size_t end_prefix( end_idx_ << (SHIFT - CMapPrefixIterator::MASK_BITS) );

    for( ; !map_prefixes.End(); map_prefixes.Next() ) {
        for( ; curr_prefix < map_prefixes.StartPrefix(); ++curr_prefix ) {
            Uint8 tval( 0 ); 
            map_file_.Write( (const char *)&tval, sizeof( Uint8 ) );
        }
        
        size_t sz( Estimate( 
                    map_prefixes.StartEntry(), map_prefixes.EndEntry() ) );

        {
            SRPRISM_ASSERT( sz < 0x100000000ULL );
            Uint8 tval( idx_file_.BytesWritten() );
            map_file_.Write( (const char *)&tval, sizeof( Uint8 ) );
            Uint4 v( (Uint4)sz );
            idx_file_.Write( (const char *)&v, sizeof( Uint4 ) );
            DBG_TRACE( "IDXOUT: bytes left: " << sz );
        }

        {
#ifndef NDEBUG
            size_t w1( idx_file_.BytesWritten() );
#endif
            FlushMapPrefix( 
                    map_prefixes.StartEntry(), map_prefixes.EndEntry() );

#ifndef NDEBUG
            SRPRISM_ASSERT( idx_file_.BytesWritten() - w1 == sz );
#endif

            ++curr_prefix;
        }

        for( ; curr_prefix < map_prefixes.EndPrefix(); ++curr_prefix ) {
            Uint8 tval( 0 ); 
            map_file_.Write( (const char *)&tval, sizeof( Uint8 ) );
        }
    }

    if( end_idx_ != NUM_HASH_KEYS ) {
        for( ; curr_prefix < end_prefix; ++curr_prefix ) {
            Uint8 tval( 0 ); 
            map_file_.Write( (const char *)&tval, sizeof( Uint8 ) );
        }
    }

    if( end_idx_ == NUM_HASH_KEYS ) {
        static const size_t NUM_MAP_PREFIX_BITS = 
            CMapPrefixIterator::MAP_PREFIX_SIZE*
            SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t NUM_MAP_PREFIXES = 
            SBitFieldTraits< size_t, NUM_MAP_PREFIX_BITS >::MAX + 1;

        for( ; curr_prefix < NUM_MAP_PREFIXES; ++curr_prefix ) {
            Uint8 tval( 0xFFFFFFFFFFFFFFFFULL );
            map_file_.Write( (const char *)&tval, sizeof( Uint8 ) );
        }
    }

    std::cout << std::flush;
}

END_NS( srprism )
END_STD_SCOPES

