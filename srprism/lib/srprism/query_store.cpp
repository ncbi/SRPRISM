/*  $Id: query_store.cpp 431273 2014-04-02 17:10:44Z morgulis $
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

#include <ncbi_pch.hpp>

#include "../common/binfile.hpp"
#include "query_store.hpp"
#include "query_acct.hpp"

START_STD_SCOPES
START_NS( srprism )

USE_NS( seq )
USE_NS( common )
    
//------------------------------------------------------------------------------
const char * CQueryStore::QDUMP_NAME = "qdump";
const char * CQueryStore::INPUT_DUMP_NAME = "batch-input-";

//------------------------------------------------------------------------------
void CQueryStore::ComputeQuerySpaceParams( 
        TSeqSize sz, int n_err, int & n_hashes, int & seed_n_err )
{
    if( sz < MIN_MED_QUERY_LEN ) { n_hashes = 2; seed_n_err = 1; }
    else if( sz < MIN_LONG_QUERY_LEN ) { n_hashes = 3; seed_n_err = 2; }
    else {
        n_hashes = std::min( MAX_N_HASHES, (int)(sz/HASH_LEN) );
        seed_n_err = std::min( (int)MAX_SEED_N_ERR, n_hashes - 1 );
    }

    seed_n_err = std::min( seed_n_err, n_err );
    n_hashes = seed_n_err + 1;
}

//------------------------------------------------------------------------------
CQueryStore::CQueryStore( 
        CMemoryManager & mem_mgr, size_t res_limit,
        TSeqSize pair_distance, TSeqSize pair_fuzz, 
        Sint2 sa_start, Sint2 sa_end, int n_err,
        bool use_fixed_hc, TWord fixed_hc )
    : mem_mgr_( mem_mgr ),
      pair_distance_( pair_distance ), pair_fuzz_( pair_fuzz ),
      sa_start_( sa_start ), sa_end_( sa_end ), n_err_( n_err ),
      size_( 0 ), res_limit_( res_limit ), n_dup_( 0 ), max_query_len_( 0 ),
      info_start_(       0 ), info_end_(       0 ),
      dup_data_start_(   0 ), dup_data_end_(   0 ),
      query_data_start_( 0 ), query_data_end_( 0 ),
      raw_data_start_(   0 ), raw_data_end_(   0 ),
      state_( INIT ),
      use_fixed_hc_( use_fixed_hc ), fixed_hc_( fixed_hc ),
      scoring_data_( 0 )
{
    if( sa_start > 0 && sa_end > 0 ) { 
        sa_start_ = sa_start - 1; 
        sa_end_ = sa_end; 
        qrv_ = false;
    }
    else if( sa_start < 0 && sa_end < 0 ) {
        sa_start_ = -sa_start - 1;
        sa_end_ = -sa_end;
        qrv_ = true;
    }
    else SRPRISM_ASSERT( false );

    SRPRISM_ASSERT( pair_fuzz_ <= pair_distance_ );
}

//------------------------------------------------------------------------------
void CQueryStore::SaveQueryData( common::CTmpStore & tmpstore )
{
    typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
    static const size_t WLETTERS = sizeof( TWord )*TTraits::PACK_FACTOR;
    static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;

    CWriteBinFile qdump( tmpstore.Register( QDUMP_NAME ) );
    query_data_end_ = query_data_start_ + size_;
    size_t n_written( 0 ), bytes_written( 0 );

    for( CQueryData * qdi( query_data_start_ ); 
            qdi != query_data_end_ && !qdi->Ignored(); ) {
        CQueryData * qdie( qdi );
        for( ; qdie != query_data_end_ && 
               CQueryData::CEqualRaw()( *qdi, *qdie ); ++qdie );
        SRPRISM_ASSERT( qdie - qdi >= 1 );
        qdump.Write( (const char *)qdi, sizeof( CQueryData ) );
        bytes_written += sizeof( CQueryData );
        size_t n_words( 3 + ((qdi->Len() - 1)>>WSHIFT) );
        n_words <<= 1;
        qdump.Write( 
                (const char *)(qdi->RawData()), n_words*sizeof( TWord ) );
        bytes_written += n_words*sizeof( TWord );
        qdi->SetPrefix( 0 ); ++qdi;
        for( ; qdi != qdie; ++qdi ) qdi->SetPrefix( 1 );
        ++n_written;
    }

    M_TRACE( CTracer::INFO_LVL, 
             "query data saved: " << n_written << " unique queries; " <<
             bytes_written << " bytes" );
    size_t qdsz( query_data_end_ - query_data_start_ );
    query_data_start_ = (CQueryData *)mem_mgr_.Shrink( 
            query_data_start_, qdsz*sizeof( CQueryData ) );
    query_data_end_ = query_data_start_ + qdsz;
}

//------------------------------------------------------------------------------
size_t CQueryStore::GenerateDuplicateData( void )
{
    size_t res( 0 );
    size_t free_space( mem_mgr_.GetFreeSpace() );
    size_t tail( free_space%sizeof( TQNum ) );
    dup_data_start_ = (TQNum *)mem_mgr_.Allocate( free_space );
    dup_data_end_ = dup_data_start_;

    for( CQueryData * qdi( query_data_start_ ); 
            qdi != query_data_end_ && !qdi->Ignored(); ) {
        CQueryData * qdie( qdi + 1 );
        while( qdie != query_data_end_ && qdie->Prefix() == 1 ) ++qdie;

        SRPRISM_ASSERT( qdie - qdi >= 1 );

        if( qdie - qdi > 1 ) {
            SRPRISM_ASSERT( 
                    free_space >= 
                        ((qdie - qdi) + DUP_HEADER_SIZE)*sizeof( TQNum ) );
            dup_data_end_[QDATA] = qdi->QNum();
            dup_data_end_[N_NOT_DONE] = dup_data_end_[N_DUP] = 
                (TQNum)(qdie - qdi);
            dup_data_end_[N_MIN_RANK] = 0;
            dup_data_end_[MIN_RANK] = res;

            for( size_t i = 0; i < dup_data_end_[N_DUP]; ++i ) {
                dup_data_end_[DUP_DATA_START + i] = qdi[i].QNum();
            }

            dup_data_end_ += (qdie - qdi) + DUP_HEADER_SIZE;
            free_space -= ((qdie - qdi) + DUP_HEADER_SIZE)*sizeof( TQNum );
            ++res;
        }

        qdi = qdie;
    }

    {
        TSize sz( dup_data_end_ - dup_data_start_ );
        dup_data_start_ = (TQNum *)mem_mgr_.Shrink( 
                dup_data_start_, tail + sz*sizeof( TQNum ) );
        dup_data_end_ = dup_data_start_ + sz;
    }
    M_TRACE( CTracer::INFO_LVL, "duplicates table created" );
    return res;
}

//------------------------------------------------------------------------------
void CQueryStore::GenerateQueryInfo( void )
{
    info_start_ = (CEntry *)mem_mgr_.Allocate( size_*sizeof( CEntry ) );
    info_end_ = info_start_ + size_;
    for( CEntry * i = info_start_; i != info_end_; ++i ) *i = CEntry();

    for( CQueryData * qdi( query_data_start_ ); 
            qdi != query_data_end_ && !qdi->Ignored(); ++qdi ) {
        CEntry & e( info_start_[qdi->QNum()] );
        e.SetIgnored( false );
        e.SetLen( qdi->Len() );
    }

    mem_mgr_.Free( query_data_start_ );
    query_data_start_ = query_data_end_ = 0;

    for( TQNum * i = dup_data_start_; i != dup_data_end_; ++i ) {
        TQNum * j( i++ );
        size_t n_dup( (size_t)*i ); i = j + DUP_HEADER_SIZE - 1;

        while( n_dup-- > 0 ) {
            CEntry & e( info_start_[*++i] );
            e.SetUnique( false );
            e.SetQueryData( j - dup_data_start_ );
        }
    }

    M_TRACE( CTracer::INFO_LVL, "initial query info setup complete" );
}

//------------------------------------------------------------------------------
void CQueryStore::ReadPrimaryQueryData( CTmpStore & tmpstore )
{
    size_t free_space( mem_mgr_.GetFreeSpace() );
    query_data_start_ = query_data_end_ = 
        (CQueryData *)mem_mgr_.Allocate( free_space );
    raw_data_start_ = raw_data_end_ = 
        (TWord *)query_data_start_ + free_space/sizeof( TWord );
    --raw_data_start_; *raw_data_start_ = 0;
    CReadBinFile qdump( tmpstore.Register( QDUMP_NAME ) );

    static const size_t WLETTERS = 
        sizeof( TWord )*SCodingTraits< SEQDATA_CODING >::PACK_FACTOR;

    while( !qdump.Eof() ) {
        // for each query read its query data structure from the file
        // and then generate additional up to two query data objects for 
        // hashes 1 and 2
        //
        if( qdump.Read( (char *)query_data_end_, sizeof( CQueryData ), true ) 
                == 0 ) {
            break;
        }

        int max_n_hashes, max_seed_n_err;
        ComputeQuerySpaceParams( 
                std::min( query_data_end_->Len(), sa_end_ - sa_start_ ), 
                n_err_, max_n_hashes, max_seed_n_err );
        SRPRISM_ASSERT( free_space >= max_n_hashes*sizeof( CQueryData ) );
        free_space -= max_n_hashes*sizeof( CQueryData );
        size_t n_words = 3 + (query_data_end_->Len() - 1)/WLETTERS;
        n_words <<= 1;
        SRPRISM_ASSERT( free_space >= sizeof( TWord )*n_words );
        raw_data_start_ -= n_words;

        if( qdump.Read( 
                    (char *)raw_data_start_, 
                    sizeof( TWord )*n_words, true ) == 0 ) {
            M_THROW( CException, IO,
                     "unexpected end of file when reading from " <<
                     tmpstore.Register( QDUMP_NAME ) );
        }

        free_space -= sizeof( TWord )*n_words;
        query_data_end_->SetRawData( raw_data_start_ );
        SetQueryData( query_data_end_ );

        {
            int max_hash( query_data_end_->GetNHashes() - 1 );

            // add all unmasked hashes
            //
            for( int i( 0 ); i < max_hash; ++i ) {
                if( !query_data_end_->HashAmbig( i ) ) {
                    query_data_end_->SetCurrHashIdx( i );
                    *(query_data_end_ + 1) = *query_data_end_;
                    ++query_data_end_;
                }
            }

            if( query_data_end_->HashAmbig( max_hash ) ) {
                query_data_end_->SetIgnored( true );
            }

            query_data_end_->SetCurrHashIdx( max_hash );
            ++query_data_end_;
        }
    }

    for( CQueryData * i( query_data_start_ ); i != query_data_end_; ++i ) {
        if( i->GetSALen() >= MIN_LONG_QUERY_LEN ) {
            info_start_[i->QNum()].SetLongBit();
        }
    }

    M_TRACE( CTracer::INFO_LVL,
             "reading of primary query data complete; " <<
             query_data_end_ - query_data_start_ << " primary queries" );
}

//------------------------------------------------------------------------------
void CQueryStore::SetupCrossLinks( void )
{
    SRPRISM_ASSERT( state_ == INIT || state_ == UPDATE );
    ClearMarks();

    for( CQueryData * i = query_data_start_; i != query_data_end_; ++i ) {
        CEntry & e( info_start_[i->QNum()] );
        TQNum offset( (TQNum)(i - query_data_start_) );

        if( e.IsUnique() ) e.SetQueryData( offset );
        else {
            e.SetMark( true );
            dup_data_start_[e.QueryData()] = offset;
        }
    }
}

//------------------------------------------------------------------------------
void CQueryStore::ReadBackQueryData( CTmpStore & tmpstore )
{
    ClearMarks();
    StartUpdate();

    if( query_data_start_ == 0 ) {
        size_t free_space( mem_mgr_.GetFreeSpace() );
        query_data_start_ = query_data_end_ = 
            (CQueryData *)mem_mgr_.Allocate( free_space );
        raw_data_start_ = raw_data_end_ = 
            (TWord *)query_data_start_ + free_space/sizeof( TWord );
    }
    else{
        query_data_end_ = query_data_start_;
        raw_data_start_ = raw_data_end_;
    }
    --raw_data_start_; *raw_data_start_ = 0;
    CReadBinFile qdump( tmpstore.Register( QDUMP_NAME ) );

    static const size_t WLETTERS = 
        sizeof( TWord )*SCodingTraits< SEQDATA_CODING >::PACK_FACTOR;

    while( !qdump.Eof() ) {
        if( qdump.Read( (char *)query_data_end_, sizeof( CQueryData ), true ) 
                == 0 ) {
            break;
        }

        TQNum qnum( query_data_end_->QNum() );

        if( !IsUnique( qnum ) ) {
            dup_data_start_[info_start_[qnum].QueryData()] = 
                query_data_end_ - query_data_start_;
            SetMark( qnum, true );
        }

        size_t n_words = 3 + (query_data_end_->Len() - 1)/WLETTERS;
        n_words <<= 1;
        raw_data_start_ -= n_words;

        if( qdump.Read( 
                    (char *)raw_data_start_, 
                    sizeof( TWord )*n_words, true ) == 0 ) {
            M_THROW( CException, IO,
                     "unexpected end of file when reading from " <<
                     tmpstore.Register( QDUMP_NAME ) );
        }

        query_data_end_->SetRawData( raw_data_start_ );
        ++query_data_end_;
    }

    StartProcess();
}

//------------------------------------------------------------------------------
size_t CQueryStore::GetRepCount( const CQueryData & q, const CRMap & rmap )
{
    size_t res( 0 );
    const TWord * raw_data( q.Data() );
    int n_hashes( q.GetNHashes() );

    for( int i( 0 ); i < n_hashes; ++i ) {
        TSeqSize hash_off( q.GetHashOffset( i ) );
        TWord word( seq::GetWord< SEQDATA_CODING >( raw_data, hash_off ) ),
              rword( 0 );
        ReverseComplement< SEQDATA_CODING >( rword, word );
        word = std::min( word, rword );
        res += rmap.RepeatRank( word );
    }

    return res;
}

//------------------------------------------------------------------------------
void CQueryStore::MarkDoneSlaves( void )
{ 
    for( size_t i(0); i < size(); ++i ) {
        if( IsSlave( i ) ) SetDone4Search< false >( i ); 
    }
}

END_NS( srprism )
END_STD_SCOPES

