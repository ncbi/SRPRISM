/*  $Id: query_store.hpp 639115 2021-10-13 15:24:22Z morgulis $
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

#ifndef __SRPRISM_QUERY_STORE_HPP__
#define __SRPRISM_QUERY_STORE_HPP__

#include <cassert>

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <common/tmpstore.hpp>
#include <common/trace.hpp>
#include <seq/seqinput.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/query_data.hpp>
#include <srprism/result.hpp>
#include <srprism/rmap.hpp>
#include <srprism/query_acct.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/tmpstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/trace.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_data.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/rmap.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_acct.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CQueryStore
{
    private:

        typedef common::Uint4 Uint4;

        typedef int TState;

        static const TState INIT         = 0;
        static const TState UPDATE       = 1;
        static const TState PROCESS      = 2;
        static const TState POST_PROCESS = 3;

        // offsets in dup table header
        static const size_t QDATA          = 0;
        static const size_t N_DUP          = 1;
        static const size_t N_NOT_DONE     = 2;
        static const size_t N_MIN_RANK     = 3;
        static const size_t MIN_RANK       = 4;
        static const size_t DUP_DATA_START = 5;
        static const size_t DUP_HEADER_SIZE = DUP_DATA_START;

        class CEntry
        {
            private:

                static const size_t LONG_Q_BIT      = 6;
                static const size_t UPRES_BIT       = 5;
                static const size_t SLAVE_BIT       = 4;
                static const size_t IGNORE_BIT      = 3;
                static const size_t UNIQ_BIT        = 2;
                static const size_t MARK_BIT        = 1;
                static const size_t DONE4SEARCH_BIT = 0;

                typedef common::Uint1 TFlagsData;
                typedef CQueryData::TQLen TQLen;
                typedef common::Uint2 THashMask;

            public:

                typedef common::Uint1 TLevelCounts;

                CEntry() 
                    : query_data_( 0 ), len_( 0 ), hash_mask_( 0 ),
                      rep_hash_count_( 0 ), rep_bu_hash_( 0 ), compl_hash_count_( 0 ),
                      flags_( (((TFlagsData)1)<<IGNORE_BIT) +
                              (((TFlagsData)1)<<UNIQ_BIT) )
                {}
                
                bool IsUPRes( void ) const {
                    return common::GetBit< UPRES_BIT >( flags_ );
                }

                void SetUPRes( void ) {
                    common::AssignBit< UPRES_BIT >( flags_, true );
                }

                bool IsLong( void ) const {
                    return common::GetBit< LONG_Q_BIT >( flags_ );
                }

                void SetLongBit( void ) {
                    common::AssignBit< LONG_Q_BIT >( flags_, true );
                }

                bool IsSlave( void ) const
                { return common::GetBit< SLAVE_BIT >( flags_ ); }

                void SetSlave( bool is_slave )
                { common::AssignBit< SLAVE_BIT >( flags_, is_slave ); }

                bool IsIgnored( void ) const
                { return common::GetBit< IGNORE_BIT >( flags_ ); }

                void SetIgnored( bool ignored )
                { common::AssignBit< IGNORE_BIT >( flags_, ignored ); }

                bool IsUnique( void ) const
                { return common::GetBit< UNIQ_BIT >( flags_ ); }

                void SetUnique( bool is_unique )
                { common::AssignBit< UNIQ_BIT >( flags_, is_unique ); }

                bool IsMarked( void ) const
                { return common::GetBit< MARK_BIT >( flags_ ); }

                void SetMark( bool marked )
                { return common::AssignBit< MARK_BIT >( flags_, marked ); }

                bool Done4Search( void ) const
                { return common::GetBit< DONE4SEARCH_BIT >( flags_ ); }

                void SetDone4Search( bool done )
                { common::AssignBit< DONE4SEARCH_BIT >( flags_, done ); }

                TQNum QueryData( void ) const { return query_data_; }

                void SetQueryData( TQNum offset ) { query_data_ = offset; }

                TSeqSize Len( void ) const { return (TSeqSize)len_; }

                void SetLen( TSeqSize len )
                {
                    SRPRISM_ASSERT( (len <= common::SIntTraits< TQLen >::MAX) );
                    len_ = len;
                }

                bool HashDone( int hash_idx ) const
                { return common::GetBit( hash_idx, hash_mask_ ); }

                void SetHashDone( int hash_idx, bool done = true )
                { common::AssignBit( hash_idx, hash_mask_, done ); }

                void ClearHashUseInfo( void )
                { compl_hash_count_ = NCompleteHashes(); hash_mask_ = 0; }

                size_t NCompleteHashes( void ) const
                { return common::CountBits( hash_mask_ ); }

                THashMask GetHashMask() const { return hash_mask_; }
                void SetHashMask( THashMask hm ) { hash_mask_ = hm; }

                bool HasRepHashes() const { 
                    return rep_hash_count_ != 0 || rep_bu_hash_ != 0;
                }

                void SetRepBUHash() { rep_bu_hash_ = 1; }
                common::Uint1 GetRepBUHash() const { return rep_bu_hash_; }

                common::Uint1 GetComplHashCount() const 
                { return compl_hash_count_; }

                void IncrRepHashCount() { ++rep_hash_count_; }

                common::Uint1 GetRepHashCount() const { 
                    return rep_hash_count_; 
                }

                void SetRepHashCount( common::Uint1 hc ) { 
                    rep_hash_count_ = hc; 
                }

            private:

                TQNum query_data_;
                TQLen len_;
                THashMask hash_mask_;
                common::Uint1 rep_hash_count_;
                common::Uint1 rep_bu_hash_;
                common::Uint1 compl_hash_count_;
                TFlagsData flags_;
        };

    public:

        static const char * QDUMP_NAME;
        static const char * INPUT_DUMP_NAME;

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode MEMORY = 0;
            static const TErrorCode IO     = 1;
            static const TErrorCode FORMAT = 2;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == MEMORY )      return "low memory";
                else if( code == IO )     return "i/o error";
                else if( code == FORMAT ) return "format error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        class CEquivIterator
        {
            public:

                CEquivIterator( const CQueryStore & qs )
                    : qs_( qs ), start_( 0 ),
                      end_( qs_.query_data_start_ )
                {
                }

                bool Next( void )
                {
                    if( end_ == qs_.query_data_end_ || end_->Ignored() ) {
                        return false;
                    }

                    start_ = end_;
                    for( ; end_ != qs_.query_data_end_ && 
                           !end_->Ignored() &&
                           CQueryData::Equiv( *start_, *end_ ); ++end_ );
                    return true;
                }

                CQueryData & Start( void ) const { return *start_; }
                CQueryData & End( void ) const { return *end_; }
                size_t Size( void ) const { return end_ - start_; }

                size_t Total( void ) const 
                { return qs_.query_data_ignored_start_ - qs_.query_data_start_; }

            private:

                const CQueryStore & qs_;
                CQueryData * start_, * end_;
        };

        CQueryStore( 
                CMemoryManager & mem_mgr, size_t res_limit,
                TSeqSize pair_distance, TSeqSize pair_fuzz,
                common::Sint2 sa_start, common::Sint2 sa_end, int n_err,
                bool use_fixed_hc, TWord fixed_hc );

        ~CQueryStore()
        {
            FreeQueryData();
            delete scoring_data_;
            if( info_start_ != 0 ) mem_mgr_.Free( (void *)info_start_ );
            if( dup_data_start_ != 0 ) mem_mgr_.Free( (void *)dup_data_start_ );
            M_TRACE( common::CTracer::INFO_LVL, "query store cleanup complete" );
        }

        template< typename t_scoring_sys >
        void Init( 
                common::CTmpStore & tmpstore, const CRMap & rmap,
                seq::CSeqInput & in, size_t max_queries, size_t ambig_limit,
                Uint4 batch_oid );

    private:

        template< bool paired, typename t_scoring_sys >
        void InitPriv(
                common::CTmpStore & tmpstore,
                const CRMap & rmap, seq::CSeqInput & in, size_t max_queries,
                size_t ambig_limit, Uint4 batch_oid );

        TQNum * DupDataStart( TQNum qn ) const
        {
            SRPRISM_ASSERT( !IsUnique( qn ) );
            return dup_data_start_ + info_start_[qn].QueryData();
        }

    public:

        size_t size(void) const { return size_; }

        bool IsUnique( TQNum qn ) const { return info_start_[qn].IsUnique(); }
        bool IsIgnored( TQNum qn ) const { return info_start_[qn].IsIgnored(); }

        TQNum * DupStart( TQNum qn ) const 
        { return DupDataStart( qn ) + DUP_DATA_START; }

        TQNum * DupEnd( TQNum qn ) const
        {
            TQNum * dup_data_start( DupDataStart( qn ) );
            return dup_data_start + DUP_DATA_START + dup_data_start[N_DUP];
        }

        CQueryData & PrimaryData( TQNum qn ) const
        {
            SRPRISM_ASSERT( state_ == INIT || state_ == PROCESS );

            if( IsUnique( qn ) ) {
                return query_data_start_[info_start_[qn].QueryData()];
            }
            else return query_data_start_[DupDataStart( qn )[QDATA]];
        }

        bool EquivMates( TQNum l, TQNum r ) const
        {
            if( l == r ) return true;
            TQNum lp( GetMate( l ) ), rp( GetMate( r ) );
            if( IsUnique( lp ) || IsUnique( rp ) ) return false;
            if( IsSlave( l ) != IsSlave( r ) ) return false;
            if( IsLeft( l ) != IsLeft( r ) ) return false;
            return DupDataStart( lp ) == DupDataStart( rp );
        }

        TQNum PrimaryQNum( TQNum qn ) const
        {
            SRPRISM_ASSERT( state_ == UPDATE || state_ == POST_PROCESS );

            if( !IsUnique( qn ) ) return DupDataStart( qn )[QDATA];
            else return qn;
        }

        static bool IsLeft( TQNum qn ) { return ((qn&0x1) == 0); }

        static TQNum GetMate( TQNum qn )
        { return IsLeft( qn ) ? qn + 1 : qn - 1; }

        template< typename t_scoring, bool paired >
        size_t NRes( TQNum qn ) const
        {
            SRPRISM_ASSERT( state_ == POST_PROCESS );
            const CQueryAcct< t_scoring > * qa( 
                    GetScoringData< t_scoring >() );

            if( !IsSlave( qn ) || !qa->HasPairedResult( qn ) ) {
                return qa->NRes( qn );
            }

            return 0;
        }

        bool Done4Search( TQNum qn ) const
        { return info_start_[qn].Done4Search(); }

        template< bool paired >
        bool GroupDone4Search( TQNum qn ) const
        {
            if( info_start_[qn].IsUnique() ) return Done4Search( qn );
            return DupDataStart( qn )[N_NOT_DONE] == 0;
        }

        template< bool paired >
        void SetDone4Search( TQNum qn, bool done = true )
        {
            assert( !IsSlave( qn ) );

            if( !info_start_[qn].IsUnique() && !Done4Search( qn ) ) {
                --DupDataStart( qn )[N_NOT_DONE];
            }

            info_start_[qn].SetDone4Search( done );
        }

        template< bool paired >
        void SetGroupDone4Search( TQNum qn )
        {
            if( IsUnique( qn ) ) {
                SetDone4Search< paired >( qn );
            }
            else {
                TQNum * dqs( DupStart( qn ) ), * dqe( DupStart( qn ) );

                for( TQNum * dqi( dqs ); dqi != dqe; ++dqi ) {
                    TQNum qn( *dqi );
                    SetDone4Search< paired >( qn );
                }
            }
        }

        bool IsMarked( TQNum qn ) const 
        { return info_start_[qn].IsMarked(); }

        void SetMark( TQNum qn, bool mark )
        { info_start_[qn].SetMark( mark ); }

        bool IsSlave( TQNum qn ) const
        { return info_start_[qn].IsSlave(); }

        bool IsUPRes( TQNum qn ) const { return info_start_[qn].IsUPRes(); }
        void SetUPRes( TQNum qn ) { info_start_[qn].SetUPRes(); }

        bool IsLong( TQNum qn ) const { return info_start_[qn].IsLong(); }

        void MarkDoneSlaves( void );

        TSeqSize Len( TQNum qn ) const { return info_start_[qn].Len(); }

        void ClearHashUseInfo( void )
        {
            for( CEntry * i( info_start_ ); i != info_end_; ++i ) {
                i->ClearHashUseInfo();
            }
        }

        void ClearMarks( void )
        {
            for( CEntry * i( info_start_ ); i != info_end_; ++i ) {
                i->SetMark( false );
            }
        }

        TSeqSize PairDistance(void) const { return pair_distance_; }
        TSeqSize PairFuzz(void) const { return pair_fuzz_; }

        void SetPairDistance( TSeqSize v ) { pair_distance_ = v; }
        void SetPairFuzz( TSeqSize v ) { pair_fuzz_ = v; }

        size_t ResLimit( void ) const { return res_limit_; }
        
        template< typename t_scoring >
        void SetResLimit( size_t new_limit )
        {
            CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );

            for( CEntry * i( info_start_ ); i != info_end_; ++i ) {
                qa->Update( i - info_start_, new_limit );
            }

            qa->SetResLim( new_limit );
            res_limit_ = new_limit;
        }

        template< typename func_t >
        void ForEach( func_t f )
        {
            SRPRISM_ASSERT( state_ == UPDATE );

            for( CQueryData * i = query_data_start_; 
                    i != query_data_end_; ++i ) {
                if( !i->Ignored() ) f( *i, this );
            }
        }
        
        void SortQueries( void )
        { 
            std::sort( 
                    query_data_start_, query_data_end_, 
                    CQueryData::CCompare() ); 

            for( query_data_ignored_start_ = query_data_start_;
                    query_data_ignored_start_ != query_data_end_ &&
                        !query_data_ignored_start_->Ignored();
                    ++query_data_ignored_start_ );
        }

        void SetQueryData( const CQueryData * q )
        {
            TQNum qnum( q->QNum() );

            if( IsUnique( qnum ) ) {
                info_start_[qnum].SetQueryData( q - query_data_start_ );
            }
            else {
                dup_data_start_[info_start_[qnum].QueryData()] =
                    q - query_data_start_;
            }
        }

        void SwapQueryData( CQueryData * q1, CQueryData * q2 )
        {
            std::swap( *q1, *q2 );
            SetQueryData( q1 );
            SetQueryData( q2 );
        }

        void SetupCrossLinks( void );

        void StartUpdate( void )      { state_ = UPDATE; }
        void StartProcess( void )     { state_ = PROCESS; }
        void StartPostProcess( void ) { state_ = POST_PROCESS; }

        void CleanUpQueryData( void ) 
        { 
            for( CEntry * i( info_start_ ); i != info_end_; ++i ) {
                if( i->IsMarked() && !i->IsUnique() ) {
                    TQNum qn( i - info_start_ );
                    DupDataStart( qn )[QDATA] = qn;
                }
            }

            query_data_end_ = query_data_start_; 
            raw_data_start_ = raw_data_end_;
            --raw_data_start_; *raw_data_start_ = 0;
        }

        void FreeQueryData( void )
        { 
            if( query_data_start_ != 0 ) mem_mgr_.Free( query_data_start_ ); 
            query_data_start_ = 0;
        }

        bool AddQueryData( 
                TWord * raw_data, size_t n_words,
                CQueryData * qdata, size_t n_queries, size_t extra = 0 )
        {
            if( (size_t)((char *)raw_data_start_ - (char *)query_data_end_) <
                    n_words*sizeof( TWord ) + 
                    n_queries*sizeof( CQueryData ) + extra) {
                return false;
            }

            raw_data_start_ -= n_words;
            std::copy( raw_data, raw_data + n_words, raw_data_start_ );

            if (qdata != nullptr)
            {
                std::copy(qdata, qdata + n_queries, query_data_end_);
            }

            while( n_queries-- > 0 ) {
                (*query_data_end_++).SetRawData( raw_data_start_ );
            }

            return true;
        }

        bool HashDone( TQNum qnum, int hash_idx ) const
        { return info_start_[qnum].HashDone( hash_idx ); }

        void SetHashDone( TQNum qnum, int hash_idx, bool done = true )
        { info_start_[qnum].SetHashDone( hash_idx, done ); }

        size_t NCompleteHashes( TQNum qnum ) const
        { return info_start_[qnum].NCompleteHashes(); }

        bool HasRepHashes( TQNum qnum ) const { 
            return info_start_[qnum].HasRepHashes(); 
        }

        void SetRepBUHash( TQNum qnum ) { info_start_[qnum].SetRepBUHash(); }

        void IncrRepHashCount( TQNum qnum ) { 
            info_start_[qnum].IncrRepHashCount(); 
        }

        bool HasRepBUHash( TQNum qn ) const {
            return info_start_[qn].GetRepBUHash() != 0;
        }

        common::Uint1 GetComplHashCount( TQNum qn ) const {
            return info_start_[qn].GetComplHashCount(); 
        }

    private:

        template< typename t_scoring >
        void AdjustMinRankInfo( TQNum qn )
        {
            CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );
            TQNum dup_idx( DupDataStart( qn )[MIN_RANK] );
            TQNum count( 0 );
            TQNum * dqs( DupStart( qn ) ), * dqe( DupEnd( qn ) );

            for( TQNum * dqi( dqs ); dqi != dqe; ++dqi ) {
                if( !IsSlave( *dqi ) ) {
                    count = qa->AdjustMinRankInfo( *dqi, dup_idx, count );
                }
            }

            SRPRISM_ASSERT( count > 0 );
            DupDataStart( qn )[N_MIN_RANK] = count;
        }

    public:

        void PairedSlavesDone4Search( TQNum qn )
        {
            SRPRISM_ASSERT( !IsUnique( qn ) );
            TQNum * dqs( DupStart( qn ) ), * dqe( DupEnd( qn ) );

            for( TQNum * dqi( dqs ); dqi != dqe; ++dqi ) {
                if( IsSlave( *dqi ) ) SetDone4Search< true >( *dqi );
            }
        }

        int MinErr( const CQueryData & q ) const
        {
            size_t n_complete_hashes( NCompleteHashes( q.QNum() ) );

            if( n_complete_hashes == 0 || q.IsSALong() ) {
                return n_complete_hashes;
            }

            if( q.GetSALen() < MIN_MED_QUERY_LEN ) {
                return n_complete_hashes - 1;
            }

            return std::min( n_complete_hashes, (size_t)1 );
        }

        template< typename t_scoring >
        bool AddPairedResult(
                TQNum * qs, TQNum * qe, int min_err, 
                TSeqSize align_len_1, 
                int n_err_1, int n_gap_1, int n_del_1, int n_gopen_1,
                TSeqSize align_len_2, 
                int n_err_2, int n_gap_2, int n_del_2, int n_gopen_2 )
        {
            CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );
            TQNum qn( *qs );
            bool adjust( !IsUnique( qn ) );
            TQNum dup_idx( adjust ? DupDataStart( qn )[MIN_RANK] : 0 );

            if( qa->AddPairedResult( 
                        qs, qe, adjust, dup_idx,
                        align_len_1, n_err_1, n_gap_1, n_del_1, n_gopen_1,
                        align_len_2, n_err_2, n_gap_2, n_del_2, n_gopen_2 ) ) {
                if( !IsUnique( qn ) ) {
                    if( DupDataStart( qn )[N_MIN_RANK] != 0 ) {
                        if( adjust ) {
                            TQNum & nmr( DupDataStart( qn )[N_MIN_RANK] );
                            SRPRISM_ASSERT( nmr >= (qe - qs) );
                            nmr -= qe - qs;

                            if( nmr == 0 ) {
                                AdjustMinRankInfo< t_scoring >( qn );
                            }
                        }
                    }
                    else AdjustMinRankInfo< t_scoring >( qn );
                }

                int max_err( qa->template MaxErr< true >( qn ) );

                if( max_err < min_err ||
                        (max_err == min_err && 
                            qa->template BestLevelFull< true >( qn )) ) {
                    for( TQNum * qi( qs ); qi != qe; ++qi ) {
                        SetDone4Search< true >( *qi );
                    }
                }

                return true;
            }
            else return false;
        }

        template< typename t_scoring >
        bool AddSingleResult(
                TQNum qn, int min_err, TSeqSize align_len,
                int n_err, int n_gap, int n_del, int n_gopen )
        {
            CQueryAcct< t_scoring > * qa( GetScoringData< t_scoring >() );
            bool res( qa->AddResult( 
                        qn, align_len, n_err, n_gap, n_del, n_gopen ) );

            if( res ) {
                int max_err( qa->template MaxErr< false >( qn ) );

                if( max_err < min_err ||
                        (max_err == min_err && 
                            qa->template BestLevelFull< false >( qn )) ) {
                    SetDone4Search< false >( qn );
                }

                return true;
            }
            else return false;
        }

        CEquivIterator EquivIterator( void ) const
        { return CEquivIterator( *this ); }

        void ReadPrimaryQueryData( common::CTmpStore & tmpstore );
        void ReadBackQueryData( common::CTmpStore & tmpstore );
        void SaveQueryDataForUnpairedRun( common::CTmpStore & tmpstore );

        template< typename t_scoring > size_t MarkDoneNotFound( void );

        template< typename t_scoring, bool paired > 
        void DistributeUnpairedCounts( void );

        TSeqSize MaxQueryLen( void ) const { return max_query_len_; }

        template< typename t_scoring, bool paired > int MaxErr( TQNum qn ) const
        { return GetScoringData< t_scoring >()->template MaxErr< paired >( qn ); }

        template< typename t_scoring, bool paired > int 
        GroupMaxErr( TQNum qn ) const
        {
            if( IsUnique( qn ) || DupDataStart( qn )[N_MIN_RANK] == 0 ) {
                return GetScoringData< t_scoring >()->template MaxErr< paired >( qn );
            }
            else {
                return GetScoringData< t_scoring >()->GroupMaxErr(
                        DupDataStart( qn )[MIN_RANK] );
            }
        }

        template< typename t_scoring >
        bool HasHigherRank( TQNum qn, const CResult & r ) const
        { return GetScoringData< t_scoring >()->HasHigherRank( qn, r ); }

        template< typename t_scoring >
        bool AddResult( TQNum qn, const CResult & r )
        { return GetScoringData< t_scoring >()->AddResult( qn, r ); }

        template< typename t_scoring >
        bool DelResult( TQNum qn, const CResult & r )
        { return GetScoringData< t_scoring >()->DelResult( qn, r ); }

        template< typename t_scoring >
        bool HasPairedResult( TQNum qn ) const
        { return GetScoringData< t_scoring >()->HasPairedResult( qn ); }

        template< typename t_scoring >
        void UpdateResLimit( TQNum qn, size_t limit )
        { GetScoringData< t_scoring >()->Update( qn, limit ); }

        /*
        void PrintQueryInfo( std::ostream & os, TQNum qn )
        {
            CEntry & e( info_start_[qn] );
            os << "\nQUERY " << qn << ": \n"
               << "flags: "
               << (e.IsSlave() ? "slave," : "master,")
               << (e.IsIgnored() ? "ignored," : "valid,")
               << (e.IsUnique() ? "unique," : "duplicate,")
               << (e.IsMarked() ? "marked," : "clear,")
               << (e.Done4Search() ? "done\n" : "active\n")
               << "hash mask: "
               << (int)e.HashDone( 0 )
               << (int)e.HashDone( 1 )
               << (int)e.HashDone( 2 ) << '\n'
               << "result: " << e.NRes() << '\n'
               << "length: " << e.Len() << '\n'
               << "rank: " << (int)e.Rank() << '\n'
               << "counts: ";
            
            const CEntry::TLevelCounts * c( e.GetCounts() );

            for( int i( 0 ); i < MAX_LEVELS; ++i ) {
                os << (int)c[i] << ' ';
            }
            
            os << '\n' << std::endl;
        }
        */

        bool QueriesReversed( void ) const { return qrv_; }
        TSeqSize GetSAStart( void ) const { return sa_start_; }
        TSeqSize GetSAEnd( void ) const { return sa_end_; }

        typedef CSeqStoreBase::TPos TPos;

        bool CheckTemplateConstraints(
                TPos left_pos, TPos right_pos,
                TSeqSize left_len, TSeqSize right_len ) const
        {
            auto lpos( std::min( left_pos, right_pos ) ),
                 rpos( std::max( left_pos + left_len, right_pos + right_len ) ),
                 tlen( rpos - lpos );
            return /* left_pos <= right_pos && */
                   tlen >= pair_distance_ - pair_fuzz_ &&
                   tlen <= pair_distance_ + pair_fuzz_;
        }

    private:

        CQueryStore( const CQueryStore & );
        CQueryStore & operator=( const CQueryStore & );

        template< typename t_scoring >
        void InitialRead( 
                common::CTmpStore & tmpstore, const CRMap & rmap,
                seq::CSeqInput & in,
                size_t max_queries, size_t ambig_limit,
                char * free_space_start, size_t free_space,
                Uint4 batch_oid );

        void SaveQueryData( common::CTmpStore & tmpstore );
        size_t GenerateDuplicateData( void );
        void GenerateQueryInfo( void );
        
        template< typename t_scoring > void GenerateScoringData( size_t n_dup );

        template< bool paired > void FinalizeDuplicateData( void ) {}
        template< bool paired > void MarkSlaves( const CRMap & rmap ) {}

        size_t GetRepCount( const CQueryData & q, const CRMap & rmap );

        static void ComputeQuerySpaceParams( 
                TSeqSize sz, int n_err, int & n_hashes, int & seed_n_err );

        template< typename t_scoring >
        CQueryAcct< t_scoring > * GetScoringData( void )
        { return static_cast< CQueryAcct< t_scoring > * >( scoring_data_ ); }

        template< typename t_scoring >
        const CQueryAcct< t_scoring > * GetScoringData( void ) const
        { 
            return static_cast< const CQueryAcct< t_scoring > * >( 
                    scoring_data_ ); 
        }

        CMemoryManager & mem_mgr_;

        TSeqSize pair_distance_, pair_fuzz_;
        TSeqSize sa_start_, sa_end_;
        bool qrv_;
        int n_err_;

        size_t size_;
        size_t res_limit_;
        size_t n_dup_;
        TSeqSize max_query_len_;

        CEntry * info_start_, * info_end_;
        TQNum * dup_data_start_, * dup_data_end_;
        CQueryData * query_data_start_, * query_data_end_;
        CQueryData * query_data_ignored_start_;
        TWord * raw_data_start_, * raw_data_end_;

        TState state_;

        bool use_fixed_hc_;
        TWord fixed_hc_;

        CQueryAcct_Base * scoring_data_;
};

template<> inline bool CQueryStore::GroupDone4Search< false >( TQNum qn ) const
{ return Done4Search( qn ); }

template<> inline void CQueryStore::SetDone4Search< false >( 
        TQNum qn, bool done )
{ info_start_[qn].SetDone4Search( done ); }

template<> inline void CQueryStore::SetGroupDone4Search< false >( TQNum qn )
{ SetDone4Search< false >( qn ); }

END_NS( srprism )
END_STD_SCOPES

#include "query_store_priv.hpp"

#endif

