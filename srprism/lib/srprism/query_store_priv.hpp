/*  $Id: query_store_priv.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: query info container
 *
 */

#ifndef NCBI_CPP_TK
#include <common/trace.hpp>
#else
#include <../src/internal/align_toolbox/srprism/lib/common/trace.hpp>
#endif

START_STD_SCOPES
START_NS( srprism )

USE_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
namespace {

    class CCompareQueriesByMate
    {
        public:

            CCompareQueriesByMate( const CQueryStore & query_store )
                : query_store_( query_store )
            {}

            bool operator()( TQNum l, TQNum r ) const
            {
                bool lslave( query_store_.IsSlave( l ) );
                bool rslave( query_store_.IsSlave( r ) );

                if( lslave && !rslave ) return false;
                else if( rslave && !lslave ) return true;

                TQNum lp( CQueryStore::GetMate( l ) );
                TQNum rp( CQueryStore::GetMate( r ) );

                if( query_store_.IsUnique( lp ) ) {
                    if( query_store_.IsUnique( rp ) ) return lp < rp;
                    else return true;
                }
                else if( query_store_.IsUnique( rp ) ) return false;

                bool lleft( query_store_.IsLeft( l ) );
                bool rleft( query_store_.IsLeft( r ) );

                if( lleft && !rleft ) return true;
                else if( rleft && !lleft ) return false;

                return (query_store_.DupStart( lp ) < 
                        query_store_.DupStart( rp ));
            }

        private:

            const CQueryStore & query_store_;
    };

}

template<> inline void CQueryStore::FinalizeDuplicateData< true >( void )
{
    for( TQNum * i = dup_data_start_; i != dup_data_end_; ) {
        TQNum * j( i + DUP_HEADER_SIZE + i[N_DUP] ); i += DUP_HEADER_SIZE;
        std::sort( i, j, CCompareQueriesByMate( *this ) );
        i = j;
    }

    for( TQNum * i = dup_data_start_; i != dup_data_end_; ) {
        TQNum * j( i + DUP_HEADER_SIZE + i[N_DUP] );
        TQNum & n_not_done( i[N_NOT_DONE] ); n_not_done = 0;
        
        for( TQNum * k( i + DUP_HEADER_SIZE ); k < j; ++k ) {
            if( !IsSlave( *k ) ) ++n_not_done;
        }

        i = j;
    }

    M_TRACE( CTracer::INFO_LVL, "duplicate data table sorted" );
}

//------------------------------------------------------------------------------
template<> inline void CQueryStore::MarkSlaves< true >( const CRMap & rmap )
{
    for( size_t i = 0; i < size_; i += 2 ) {
        if( !info_start_[i].IsIgnored() && !info_start_[i+1].IsIgnored() ) {
            const CQueryData & q( PrimaryData( i ) );
            const CQueryData & m( PrimaryData( i + 1 ) );

            size_t rc1( GetRepCount( q, rmap ) );
            size_t rc2( GetRepCount( m, rmap ) );

            if( rc1 < rc2 ) info_start_[i+1].SetSlave( true );
            else info_start_[i].SetSlave( true );
        }
    }

    M_TRACE( CTracer::INFO_LVL,
             "master/slave relationship between mates established" );
}

//------------------------------------------------------------------------------
template< typename t_scoring, bool paired >
inline void CQueryStore::DistributeUnpairedCounts( void )
{
    CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );

    if( paired ) {
        for( size_t i( 0 ); i < size_; i += 2 ) {
            TQNum slave( IsSlave( i ) ? i : i + 1 ), master( GetMate( slave ) );

            if( qa->HasPairedResult( master ) ) {
                qa->DupScoringData( slave, master );
            }
        }

        for( TQNum * i = dup_data_start_; i != dup_data_end_; ) {
            TQNum * ie( i + DUP_HEADER_SIZE + i[N_DUP] );
            TQNum * di( i + DUP_DATA_START ), * j( di );
            while( j != ie && !HasRepBUHash( *j ) ) ++j;

            if( j != ie ) for( ; di != ie; ++di ) {
                if( di != j ) info_start_[*di].SetRepBUHash();
            }

            i = ie;
        }
    }
    else {
        for( TQNum * i = dup_data_start_; i != dup_data_end_; ) {
            TQNum * ie( i + DUP_HEADER_SIZE + i[N_DUP] );
            TQNum * di( i + DUP_DATA_START ), * j( di );
            while( j != ie && qa->NRes( *j ) == 0 ) ++j;

            if( j != ie ) {
                for( ; di != ie; ++di ) {
                    if( di != j ) qa->DupScoringData( *di, *j );
                }
            }

            di = i + DUP_DATA_START; 
            j = di;

            while( j != ie && 
                    !info_start_[*j].HasRepHashes() && 
                    info_start_[*j].NCompleteHashes() == 0 ) {
                ++j;
            }

            if( j != ie ) {
                CEntry & ee( info_start_[*j] );

                for( ; di != ie; ++di ) {
                    if( di != j ) {
                        CEntry & e( info_start_[*di] );
                        e.SetHashMask( ee.GetHashMask() );
                        e.SetRepHashCount( ee.GetRepHashCount() );
                        if( ee.GetRepBUHash() != 0 ) e.SetRepBUHash();
                    }
                }
            }

            di = i + DUP_DATA_START;
            j = di;
            for( ; j != ie && !info_start_[*j].IsLong(); ++j );

            if( j != ie ) {
                for( ; di != ie; ++di ) {
                    if( di != j ) info_start_[*di].SetLongBit();
                }
            }

            i = ie;
        }
    }
}

//------------------------------------------------------------------------------
template< typename t_scoring >
void CQueryStore::InitialRead( 
        common::CTmpStore & tmpstore, const CRMap & rmap, CSeqInput & in,
        size_t max_queries, size_t ambig_limit,
        char * free_space_start, size_t free_space, Uint4 batch_oid )
{
    SRPRISM_ASSERT( state_ == INIT );

    typedef SCodingTraits< SEQDATA_CODING > TTraits;
    static const size_t WLETTERS = sizeof( TWord )*TTraits::PACK_FACTOR;

    SRPRISM_ASSERT( (size_t)in.NCols() <= SIntTraits< TQNum >::MAX );
    max_queries = std::min( 
            max_queries, (size_t)SIntTraits< TQNum >::MAX - in.NCols() + 1 );
    M_TRACE( common::CTracer::INFO_LVL,
             "reading up to " << max_queries << " queries" );
    bool truncate_warning = false;
    bool ambig_warning    = false;
    bool short_warning    = false;
    size_t ignored = 0, truncated = 0, ambig = 0, too_short = 0;
    CQueryData * qdata_end( (CQueryData *)free_space_start ),
               * qdata_start( qdata_end );
    TWord * qraw_start( 
            (TWord *)free_space_start + free_space/sizeof( TWord ) );
    free_space = ((char *)qraw_start - (char *)qdata_end );

    {
        int n_cols( in.NCols() );

        //######################################################################
        //
        // Conservative estimate of per-query information stored in memory.
        //
        size_t acct_bpq( CQueryAcct< t_scoring >::EstimateBytesPerQuery(
                    res_limit_, n_err_, (n_cols == 2) ) ); // size of score
                                                           // accounting data
        size_t qsz_estimate( 
                sizeof( CEntry ) + // size of query info in query store
                MAX_N_HASHES*sizeof( CQueryData ) + // size of seed data
                2*MAX_QUERY_LEN*TTraits::PACK_FACTOR + // size of sequence data
                5*sizeof( TWord ) + // sentinel data on both sides of 
                                    // sequence data and one final sentinel
                acct_bpq );
        qsz_estimate *= n_cols;
        //######################################################################

        std::string input_dump_name( INPUT_DUMP_NAME );
        input_dump_name += std::to_string( batch_oid );
        CWriteTextFile_CPPStream idump( tmpstore.Register( input_dump_name ) );

        for( size_t i = 0; i < max_queries && !in.Done(); ) {
            if( qsz_estimate > free_space ) break;
            if( !in.Next() ) break;
            free_space -= n_cols*(sizeof( CEntry ) + acct_bpq);
            TSeqId id( in.Id() );

            for( size_t j = 0; j < (size_t)n_cols; ++j, ++i ) {
                typedef CSeqInput::TData TSrcData;
                const TSrcData & data( in.Data( j ) );

                idump.LineOut( std::string( ">" ) + id );
                idump.LineOut( std::string(
                    data.seq.begin(), data.seq.begin() + data.size ) );

                bool ignore = false;
                TSeqSize data_size( data.size );
                
                if( data_size > MAX_QUERY_LEN ) {
                    if( !truncate_warning ) {
                        M_TRACE( common::CTracer::WARNING_LVL, 
                                 "sequences longer than " << MAX_QUERY_LEN << 
                                 " bases will be truncated" );
                        truncate_warning = true;
                    }

                    ++truncated;
                    data_size = MAX_QUERY_LEN;
                }

                if( data_size < HASH_LEN ) {
                    if( !short_warning ) {
                        M_TRACE( common::CTracer::WARNING_LVL,
                                 "sequences shorter than " << HASH_LEN << 
                                 " bases will be ignored" );
                        short_warning = true;
                    }

                    ignore = true;
                    ++too_short;
                }

                size_t n_ambig( data.Ambig( sa_start_, sa_end_ ) );
    
                if( n_ambig > ambig_limit ) {
                    if( !ambig_warning ) {
                        M_TRACE( common::CTracer::WARNING_LVL,
                                 "sequences with over " << ambig_limit << 
                                 " ambiguities in the seeding area will be "
                                 " ignored" );
                        ambig_warning = true;
                    }

                    ignore = true;
                    ++ambig;
                }

                {
                    int max_n_hashes, max_seed_n_err;
                    ComputeQuerySpaceParams( 
                            std::min( data_size, sa_end_ ) - sa_start_,
                            n_err_, max_n_hashes, max_seed_n_err );
                    TWord * raw_data( 0 );

                    if( data_size > 0 ) {
                        size_t n_words( 3 + ((data_size - 1)/WLETTERS));
                        n_words <<= 1;
                        qraw_start -= n_words;
                        raw_data = qraw_start;

                        TSrcData::TSeq rvsec;

                        if( qrv_ ) {
                            for( size_t j( 0 ); j < data_size; ++j ) {
                                rvsec.push_back( 
                                        SCodingTraits< TSrcData::CODING >::
                                            RC[data.seq[data_size - j - 1]] );
                            }
                        }

                        TSrcData::TSeq & s( qrv_ ? rvsec : data.seq );

                        std::fill( raw_data, raw_data + n_words, 0 );
                        Recode< SEQDATA_CODING, TSrcData::CODING >(
                                raw_data + 2, &s[0], data_size );

                        {
                            TWord * ambig_data( raw_data + (n_words>>1) );
                            std::fill( 
                                    ambig_data, ambig_data + (n_words>>1), 
                                    0xFFFFFFFF );

                            for( size_t j( 0 ); j < data_size; ++j ) {
                                if( SCodingTraits< TSrcData::CODING >::
                                        NAMBIG[s[j]] ) {
                                    SetStreamLetter< SEQDATA_CODING >(
                                            ambig_data + 1, j, 0 );
                                }
                                else {
                                    SetStreamLetter< SEQDATA_CODING >(
                                            ambig_data + 1, j, 3 );
                                }
                            }
                        }

                        free_space -= n_words*sizeof( TWord );
    
                        if( data_size > max_query_len_ ) {
                            max_query_len_ = data_size;
                        }
                    }

                    *qdata_end = CQueryData( size_++, data_size, raw_data );
                    qdata_end->SetAmbig( n_ambig > 0 );

                    if( !ignore ) {
                        ignore = !qdata_end->ComputeHC( 
                                rmap, sa_start_, sa_end_, max_seed_n_err,
                                use_fixed_hc_, fixed_hc_ );
                    }

                    if( ignore ) ++ignored;
                    qdata_end->SetIgnored( ignore );
                    ++qdata_end;
                    free_space -= 
                        max_n_hashes*sizeof( CQueryData ) + sizeof( TWord );
                }
            }
        }
    }

    SRPRISM_ASSERT( size_ == (size_t)(qdata_end - qdata_start) );
    M_TRACE( CTracer::INFO_LVL, 
             "got " << qdata_end - qdata_start << " queries" );
    std::sort( qdata_start, qdata_end, CQueryData::CCompareRaw() );
    query_data_start_ = qdata_start;
    M_TRACE( CTracer::INFO_LVL, size_ << " queries read" );
}

//------------------------------------------------------------------------------
template< typename t_scoring > 
void CQueryStore::GenerateScoringData( size_t n_dup )
{
    scoring_data_ = new CQueryAcct< t_scoring >( 
            mem_mgr_, size_, res_limit_, n_dup, n_err_ );
}

//------------------------------------------------------------------------------
template< typename t_scoring_sys >
void CQueryStore::Init( 
        CTmpStore & tmpstore, const CRMap & rmap, CSeqInput & in, 
        size_t max_queries, size_t ambig_limit, Uint4 batch_oid )
{
    switch( in.NCols() ) {
        case 1: 
            InitPriv< false, t_scoring_sys >( 
                    tmpstore, rmap, in, max_queries, ambig_limit, batch_oid );
            break;

        case 2:
            InitPriv< true, t_scoring_sys >( 
                    tmpstore, rmap, in, max_queries, ambig_limit, batch_oid );
            break;

        default:

            M_THROW( CException, FORMAT, 
                     "unexpected number of input columns: " << in.NCols() );
    }
}

//------------------------------------------------------------------------------
template< bool paired, typename t_scoring_sys >
void CQueryStore::InitPriv( 
        CTmpStore & tmpstore, const CRMap & rmap, 
        CSeqInput & in, size_t max_queries, size_t ambig_limit,
        Uint4 batch_oid )
{
    size_t free_space( mem_mgr_.GetFreeSpaceSize() );
    char * free_space_start = 
        (char *)mem_mgr_.Allocate( mem_mgr_.GetFreeSpaceSize() );

    try {
        InitialRead< t_scoring_sys >( 
                tmpstore, rmap, in, max_queries, ambig_limit,
                free_space_start, free_space, batch_oid );
    }
    catch( ... ) {
        mem_mgr_.Free( free_space_start );
        throw;
    }

    M_TRACE( CTracer::INFO_LVL, "initial query sort complete" );
    SaveQueryData( tmpstore );
    size_t n_dup( GenerateDuplicateData() );
    GenerateQueryInfo();
    GenerateScoringData< t_scoring_sys >( n_dup );
    ReadPrimaryQueryData( tmpstore );
    MarkSlaves< paired >( rmap );
    FinalizeDuplicateData< paired >();
}

//------------------------------------------------------------------------------
template< typename t_scoring >
size_t CQueryStore::MarkDoneNotFound( void )
{
    const CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );
    size_t res( 0 );

    for( size_t i(0); i < size(); i += 2 ) {
        if( qa->NRes( i ) == 0 && qa->NRes( i + 1 ) == 0 ) {
            SetDone4Search< false >( i );
            SetDone4Search< false >( i + 1 );
            ++res;
        }
    }

    return res;
}

END_NS( srprism )
END_STD_SCOPES

