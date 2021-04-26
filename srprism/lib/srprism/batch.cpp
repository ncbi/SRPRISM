/*  $Id: batch.cpp 573027 2018-10-22 14:43:30Z morgulis $
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
 * File Description: srprism query batch processor
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <set>

#include "batch.hpp"
#include "search_mode.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
const char * CBatch::TMP_RES_FNAME = "tmpres";

const double CBatch::ISD_ACCEPT_PER_SC_RATIO = 0.985;
const double CBatch::ISD_ACCEPT_RATIO = 0.9;

//------------------------------------------------------------------------------
CBatch::CBatch( 
        SBatchInitData & init_data, 
        CSeqInput & in, TQueryOrdId start_qid )
    : init_data_( init_data ),
      tmp_store_( init_data.tmpdir ),
      /*
      u_tmpres_mgr_( 
              init_data.u_tmp_res_buf, init_data.u_tmp_res_buf_size, 
              "utmpres", tmp_store_ ),
      p_tmpres_mgr_( 
              init_data.p_tmp_res_buf, init_data.p_tmp_res_buf_size, 
              "ptmpres", tmp_store_ ),
      tmpres_mgr_p_( &p_tmpres_mgr_ ),
      */
      seqstore_( *init_data.seqstore_p ), rmap_( init_data.index_basename ),
      use_sids_( init_data.use_sids ), use_qids_( init_data.use_qids ),
      search_mode_( init_data.search_mode ),
      final_res_limit_( init_data.res_limit + 1 ),
      start_qid_( start_qid ), end_qid_( start_qid ), queries_p_( 0 ),
      out_p_( init_data.out_p ),
      paired_log_( init_data.paired_log )
{
    seqstore_.Load();

    if( init_data_.n_threads > 1 )
    {
        auto free_space(
            (init_data_.mem_mgr_p->GetFreeSpace())/init_data_.n_threads );
        init_data_.mem_mgr_p.reset( new CMemoryManager( free_space ) );
        char * t( (char *)init_data_.mem_mgr_p->Allocate( TMP_RES_BUF_SIZE ) );
        init_data_.u_tmp_res_buf = t;
        t = (char *)init_data_.mem_mgr_p->Allocate( TMP_RES_BUF_SIZE );
        init_data_.p_tmp_res_buf = t;
    }

    Uint4 u_res_limit( 0 );
    
    if( search_mode_ == SSearchMode::DEFAULT ) {
        typedef SSearchModeTraits< SSearchMode::DEFAULT >::TScoringSys 
            TScoringSys;

        res_limit_ = std::max( 
                (Uint4)MIN_RES_LIMIT,
                final_res_limit_*seqstore_.OverlapFactor() );
        res_limit_ = std::min( res_limit_, (Uint4)TScoringSys::MAX_N_RES );

        u_res_limit = 
                init_data.paired 
                    ? std::min( (Uint4)TScoringSys::MAX_N_RES, 5*res_limit_ ) 
                    : res_limit_;

        queries_p_.reset( new CQueryStore( 
                    *init_data.mem_mgr_p, u_res_limit, 
                    init_data.pair_distance, init_data.pair_fuzz, 
                    init_data.sa_start, init_data.sa_end, init_data.n_err, 
                    init_data.use_fixed_hc, (TWord)init_data.fixed_hc ) );
        queries_p_->Init< TScoringSys >( 
                tmp_store_, rmap_, in, 
                (size_t)init_data.batch_limit, init_data.n_err );
        queries_p_->SetResLimit< TScoringSys >( u_res_limit );
    }
    else if( search_mode_ == SSearchMode::PARTIAL ) {
        typedef SSearchModeTraits< SSearchMode::PARTIAL >::TScoringSys 
            TScoringSys;

        res_limit_ = std::max( 
                (Uint4)MIN_RES_LIMIT,
                final_res_limit_*seqstore_.OverlapFactor() );
        res_limit_ = std::min( res_limit_, (Uint4)TScoringSys::MAX_N_RES );

        u_res_limit = 
                init_data.paired 
                    ? std::min( (Uint4)TScoringSys::MAX_N_RES, 5*res_limit_ ) 
                    : res_limit_;

        queries_p_.reset( new CQueryStore( 
                    *init_data.mem_mgr_p, u_res_limit, 
                    init_data.pair_distance, init_data.pair_fuzz, 
                    init_data.sa_start, init_data.sa_end, init_data.n_err, 
                    init_data.use_fixed_hc, (TWord)init_data.fixed_hc ) );
        queries_p_->Init< TScoringSys >( 
                tmp_store_, rmap_, in, 
                (size_t)init_data.batch_limit, init_data.n_err );
        queries_p_->SetResLimit< TScoringSys >( u_res_limit );
    }
    else if( search_mode_ == SSearchMode::SUM_ERR ) {
        typedef SSearchModeTraits< SSearchMode::SUM_ERR >::TScoringSys 
            TScoringSys;

        res_limit_ = std::max( 
                (Uint4)MIN_RES_LIMIT,
                final_res_limit_*seqstore_.OverlapFactor() );
        res_limit_ = std::min( res_limit_, (Uint4)TScoringSys::MAX_N_RES );

        u_res_limit =  init_data.paired ? std::min( 
                                            (Uint4)TScoringSys::MAX_N_RES, 
                                            5*res_limit_ ) 
                                        : res_limit_;

        queries_p_.reset( 
                new CQueryStore( 
                    *init_data.mem_mgr_p, u_res_limit,
                    init_data.pair_distance, init_data.pair_fuzz,
                    init_data.sa_start, init_data.sa_end, init_data.n_err,
                    init_data.use_fixed_hc, (TWord)init_data.fixed_hc ) );
        queries_p_->Init< TScoringSys >( 
                tmp_store_, rmap_, in, 
                (size_t)init_data.batch_limit, init_data.n_err );
        queries_p_->SetResLimit< TScoringSys >( u_res_limit );
    }
    else if( search_mode_ == SSearchMode::BOUND_ERR ) {
        typedef SSearchModeTraits< SSearchMode::BOUND_ERR >::TScoringSys 
            TScoringSys;

        res_limit_ = std::max( 
                (Uint4)MIN_RES_LIMIT,
                final_res_limit_*seqstore_.OverlapFactor() );
        res_limit_ = std::min( res_limit_, (Uint4)TScoringSys::MAX_N_RES );

        u_res_limit =  init_data.paired ? std::min( 
                                            (Uint4)TScoringSys::MAX_N_RES, 
                                            5*res_limit_ ) 
                                        : res_limit_;

        queries_p_.reset( 
                new CQueryStore( 
                    *init_data.mem_mgr_p, u_res_limit,
                    init_data.pair_distance, init_data.pair_fuzz,
                    init_data.sa_start, init_data.sa_end, init_data.n_err,
                    init_data.use_fixed_hc, (TWord)init_data.fixed_hc ) );
        queries_p_->Init< TScoringSys >( 
                tmp_store_, rmap_, in, 
                (size_t)init_data.batch_limit, init_data.n_err );
        queries_p_->SetResLimit< TScoringSys >( u_res_limit );
    }
    else SRPRISM_ASSERT( false );

    end_qid_ = start_qid + queries_p_->size();

    // prepare pass initialization data
    {
        pass_init_data_.search_stats = init_data.search_stats;
        pass_init_data_.index_basename = init_data.index_basename;
        pass_init_data_.res_limit = u_res_limit;
        pass_init_data_.repeat_threshold = init_data.repeat_threshold;
        pass_init_data_.pair_distance = init_data.pair_distance;
        pass_init_data_.pair_fuzz = init_data.pair_fuzz;
        pass_init_data_.n_err = init_data.n_err;

        pass_init_data_.ipam_vec = init_data.ipam_vec;

        // pass_init_data_.mem_mgr_p = init_data.mem_mgr_p;
        pass_init_data_.mem_mgr_p = init_data.mem_mgr_p.get();
        pass_init_data_.seqstore_p = &seqstore_;
        pass_init_data_.tmp_store_p = &tmp_store_;
        pass_init_data_.randomize = init_data.randomize;
        pass_init_data_.random_seed = init_data.random_seed;

        pass_init_data_.queries_p = queries_p_.get();
    }
}

//------------------------------------------------------------------------------
template<> bool CBatch::Run< true >( void )
{
    pass_init_data_.paired_search = true;
    // pass_init_data_.tmpres_mgr_p = &u_tmpres_mgr_;
    pass_init_data_.tmpres_mgr_p = u_tmpres_mgr_.get();

    if( pass_init_data_.n_err > 2 ) {
        RunPass< HASH_NORMAL, false >( false, true );
        RunPass< HASH_NORMAL, false >( true, false );
    }
    else RunPass< HASH_NORMAL, false >( false, false );

    RunPass< HASH_BLOWUP, false >( false, false );
    bool cont( true );

    if( search_mode_ == SSearchMode::DEFAULT ) {
        cont = InterProcess< SSearchMode::DEFAULT >();
    }
    else if( search_mode_ == SSearchMode::SUM_ERR ) {
        cont = InterProcess< SSearchMode::SUM_ERR >();
    }
    else if( search_mode_ == SSearchMode::PARTIAL ) {
        cont = InterProcess< SSearchMode::PARTIAL >();
    }
    else if( search_mode_ == SSearchMode::BOUND_ERR ) {
        InterProcess< SSearchMode::BOUND_ERR >();
    }
    else SRPRISM_ASSERT( false );

    if( !cont ) {
        // tmpres_mgr_p_ = &u_tmpres_mgr_;
        tmpres_mgr_p_ = u_tmpres_mgr_.get();

        if( search_mode_ == SSearchMode::DEFAULT ) {
            PostProcess< SSearchMode::DEFAULT, false >();
        }
        else if( search_mode_ == SSearchMode::SUM_ERR ) {
            PostProcess< SSearchMode::SUM_ERR, false >();
        }
        else if( search_mode_ == SSearchMode::PARTIAL ) {
            PostProcess< SSearchMode::PARTIAL, false >();
        }
        else SRPRISM_ASSERT( false );

        return false;
    }

    // pass_init_data_.tmpres_mgr_p = &p_tmpres_mgr_;
    pass_init_data_.tmpres_mgr_p = p_tmpres_mgr_.get();
    queries_p_->SetPairDistance( init_data_.pair_distance );
    queries_p_->SetPairFuzz( init_data_.pair_fuzz );

    if( pass_init_data_.n_err > 1 ) {
        RunPass< HASH_NORMAL, true >( false, true );
        RunPass< HASH_NORMAL, true >( true, false );
    }
    else RunPass< HASH_NORMAL, true >( false, false );

    RunPass< HASH_BLOWUP, true >( false, false );

    if( search_mode_ == SSearchMode::DEFAULT ) {
        PostProcess< SSearchMode::DEFAULT, true >();
    }
    else if( search_mode_ == SSearchMode::SUM_ERR ) {
        PostProcess< SSearchMode::SUM_ERR, true >();
    }
    else if( search_mode_ == SSearchMode::PARTIAL ) {
        PostProcess< SSearchMode::PARTIAL, true >();
    }
    else if( search_mode_ == SSearchMode::BOUND_ERR ) {
        PostProcess< SSearchMode::BOUND_ERR, true >();
    }
    else SRPRISM_ASSERT( false );

    return true;
}

//------------------------------------------------------------------------------
template<> bool CBatch::Run< false >( void )
{
    pass_init_data_.paired_search = false;
    // pass_init_data_.tmpres_mgr_p = &p_tmpres_mgr_;
    pass_init_data_.tmpres_mgr_p = p_tmpres_mgr_.get();

    if( pass_init_data_.n_err > 2 ) {
        RunPass< HASH_NORMAL, false >( false, true );
        RunPass< HASH_NORMAL, false >( true, false );
    }
    else RunPass< HASH_NORMAL, false >( false, false );

    RunPass< HASH_BLOWUP, false >( false, false );

    if( search_mode_ == SSearchMode::DEFAULT ) {
        PostProcess< SSearchMode::DEFAULT, false >();
    }
    else if( search_mode_ == SSearchMode::SUM_ERR ) {
        PostProcess< SSearchMode::SUM_ERR, false >();
    }
    else if( search_mode_ == SSearchMode::PARTIAL ) {
        PostProcess< SSearchMode::PARTIAL, false >();
    }
    else if( search_mode_ == SSearchMode::BOUND_ERR ) {
        PostProcess< SSearchMode::BOUND_ERR, false >();
    }
    else SRPRISM_ASSERT( false );

    return true;
}

//------------------------------------------------------------------------------
void CBatch::DiversifySubjects( CResult * s, CResult * e )
{
    if( e -s < 2 ) return;
    std::set< TDBOrdId > sids;
    for( CResult * i( s ); i != e; ++i ) sids.insert( i->SNum() );
    size_t j( 0 );

    for( CResult * i( s ); i != e - j; ) {
        std::set< TDBOrdId >::iterator sidi( sids.find( i->SNum() ) );

        if( sidi != sids.end() ) { sids.erase( sidi ); ++i; }
        else { ++j; std::swap( *i, *(e-j) ); }
    }
}

//------------------------------------------------------------------------------
namespace {
    struct SResultCompare
    {
        bool operator()( const CResult & l, const CResult & r ) const
        {
            if( l.SNum() == r.SNum() ) return l.SOff( 0 ) < r.SOff( 0 );
            else return l.SNum() < r.SNum();
        }
    };
}

//------------------------------------------------------------------------------
inline int CBatch::GetStrandConfig( CResult * l, CResult * r )
{
    if( l->Strand( 0 ) == STRAND_FW ) {
        if( l->SOff( 0 ) < r->SOff( 0 ) ) {
            return (r->Strand( 0 ) == STRAND_FW) ? 0 : 1;
        }
        else return (r->Strand( 0 ) == STRAND_FW ) ? 2 : 3;
    }
    else {
        if( l->SOff( 0 ) < r->SOff( 0 ) ) {
            return (r->Strand( 0 ) == STRAND_FW) ? 3 : 2;
        }
        else return (r->Strand( 0 ) == STRAND_FW) ? 1 : 0;
    }
}

//------------------------------------------------------------------------------
inline bool CBatch::CheckStrandConfig( CResult * l, CResult * r )
{
    size_t ipam_idx( l->Strand( 0 ) == STRAND_FW ? 0 : 1 );
    const T_IPAM & ipam( pass_init_data_.ipam_vec.data[ipam_idx] );
    // T_IPAM lcl_ipam( 0xF );
    T_IPAM lcl_ipam( 0 );

    /*
    if( l->SOff( 0 ) < r->SOff( 0 ) ) lcl_ipam = (ipam&IPAM_RIGHT_ENABLED);
    else if( l->SOff( 0 ) > r->SOff( 0 ) ) lcl_ipam = (ipam&IPAM_LEFT_ENABLED);
    else return false;
    */
    if( l->SOff( 0 ) <= r->SOff( 0 ) ) lcl_ipam |= (ipam&IPAM_RIGHT_ENABLED);
    if( l->SOff( 0 ) >= r->SOff( 0 ) ) lcl_ipam |= (ipam&IPAM_LEFT_ENABLED);

    if( r->Strand( 0 ) == STRAND_FW ) return (lcl_ipam&IPAM_FW_ENABLED) != 0;
    else return (lcl_ipam&IPAM_RV_ENABLED) != 0;
}

//------------------------------------------------------------------------------
inline void CBatch::CheckPair( 
        CResult * l, CResult * r, TSeqSize e, size_t & res,
        TPairCandidates & pc )
{
    TSeqSize e1( r->SOff( 0 ) + r->GetAlignLen( 0 ) - r->GetNIns( 0 ) );

    if( e1 <= e && CheckStrandConfig( l, r ) ) {
        pc.push_back( std::make_pair( l, r ) );
        ++res;
    }
}

//------------------------------------------------------------------------------
size_t CBatch::IdentifyPairsForSubject(
        CResult * ls, CResult * le,
        CResult * rs, CResult * re,
        TPairCandidates & pc )
{
    size_t res( 0 );
    TSeqSize d( init_data_.pair_distance ),
             f( init_data_.pair_fuzz );
    SRPRISM_ASSERT( f <= d );

    for( ; ls < le && rs < re; ++ls ) {
        TSeqSize lsoff( ls->SOff(0) );
        TSeqSize ee( lsoff + ls->GetAlignLen( 0 ) - ls->GetNIns( 0 ) ), 
                 e( ee - std::min( ee, d - f ) ),
                 s( ee - std::min( ee, d + f ) );
        while( rs < re && rs->SOff( 0 ) < s ) ++rs;
        CResult * r( rs );

        // check to the left of ls
        //
        for( ; r < re && r->SOff( 0 ) < e; ++r ) {
            CheckPair( ls, r, ee, res, pc );
        }

        s = lsoff + (d - f); e = lsoff + d + f; r = rs;
        while( r < re && r->SOff( 0 ) < lsoff ) ++r;

        // check to the right of ls
        //
        for( ; r < re && r->SOff( 0 ) < e; ++r ) {
            if( r->SOff( 0 ) + r->GetAlignLen( 0 ) - r->GetNIns( 0 ) >= s ) {
                CheckPair( ls, r, e, res, pc );
            }
        }
    }

    return res;
}

//------------------------------------------------------------------------------
template< int search_mode >
size_t CBatch::IdentifyPairs( CResult * s, size_t n_left, size_t n_right )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    TQNum qn( s->QNum() );
    if( queries_p_->IsSlave( qn ) ) qn = queries_p_->GetMate( qn );
    TPairCandidates pair_candidates;
    CResult * ls( s ),  * le( ls + n_left ),
            * rs( le ), * re( rs + n_right ),
            * lss( ls ), * rss( rs ), * lee( lss ), * ree( rss );
    std::stable_sort( ls, le, SResultCompare() );
    std::stable_sort( rs, re, SResultCompare() );
    size_t res( 0 );

    while( lss < le && rss < re ) {
        TDBOrdId lsnum( lss->SNum() );
        while( lee != le && lee->SNum() == lsnum ) ++lee;
        while( rss != re && rss->SNum() < lsnum ) { ++rss; ++ree; }

        if( rss != re && rss->SNum() == lsnum ) {
            while( ree != re && ree->SNum() == lsnum ) ++ree;
            res += IdentifyPairsForSubject( 
                    lss, lee, rss, ree, pair_candidates );
            rss = ree;
        }

        lss = lee;
    }

    for( TPairCandidates::const_iterator i( pair_candidates.begin() );
            i != pair_candidates.end(); ++i ) {
        if( !seqstore_.CheckRegionPair(
                    seqstore_.EncodePos( 
                        std::make_pair( 
                            i->first->SNum(), 
                            i->first->SOff( 0 ) + i->first->GetLeftOffset( 0 ) ) ),
                    seqstore_.EncodePos( 
                        std::make_pair( 
                            i->second->SNum(), 
                            i->second->SOff( 0 ) + i->second->GetLeftOffset( 0 ) ) ),
                    i->first->GetAlignLen( 0 ) - i->first->GetNIns( 0 ),
                    i->second->GetAlignLen( 0 ) - i->second->GetNIns( 0 ) ) ) {
            continue;
        }

        if( queries_p_->AddPairedResult< TScoring >( 
                    &qn, &qn + 1, 0,
                    i->first->GetAlignLen( 0 ),
                    i->first->NErr( 0 ), 
                    i->first->GetNId( 0 ), 
                    i->first->GetNDel( 0 ),
                    i->first->GetNGOpen( 0 ),
                    i->second->GetAlignLen( 0 ),
                    i->second->NErr( 0 ), 
                    i->second->GetNId( 0 ), 
                    i->second->GetNDel( 0 ),
                    i->second->GetNGOpen( 0 ) ) ) {
            // CResult r( p_tmpres_mgr_.Save( CResult::EstimateLen( 
            CResult r( p_tmpres_mgr_->Save( CResult::EstimateLen( 
                            2, i->first->NErr( 0 ), i->second->NErr( 0 ) ) ) );
            r.Init( 
                    qn, i->first->SNum(), 
                    i->first->RawSOff( 0 ), i->second->RawSOff( 0 ),
                    i->first->GetAlignLen( 0 ), i->second->GetAlignLen( 0 ),
                    i->first->GetRawLeftOffset( 0 ), 
                    i->second->GetRawLeftOffset( 0 ),
                    i->first->GetRawRightOffset( 0 ), 
                    i->second->GetRawRightOffset( 0 ),
                    i->first->Strand( 0 ), i->second->Strand( 0 ),
                    i->first->NErr( 0 ), i->second->NErr( 0 ) );

            if( i->first->NErr( 0 ) > 0 ) {
                CResult::CErrorIterator e( i->first->ErrorIterator( 0 ) );
                size_t c( 0 );

                do {
                    r.SetErrorInfo( 0, c++, e.QOff(), e.ErrType() );
                    e.Next();
                } while( !e.End() );
            }

            if( i->second->NErr( 0 ) > 0 ) {
                CResult::CErrorIterator e( i->second->ErrorIterator( 0 ) );
                size_t c( 0 );

                do {
                    r.SetErrorInfo( 1, c++, e.QOff(), e.ErrType() );
                    e.Next();
                } while( !e.End() );
            }
        }
    }

    return res;
}

//------------------------------------------------------------------------------
template< int search_mode >
void CBatch::UpdateCounts( CResult * s, CResult * e )
{ 
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;

    for( ; s < e; ++s ) {
        queries_p_->UpdateResLimit< TScoring >( s->QNum(), res_limit_ );
    }
}

//------------------------------------------------------------------------------
void CBatch::MarkDone( CResult * s, CResult * e )
{ 
    for( ; s < e ; ++s ) {
        TQNum qn( s->QNum() );

        if( !queries_p_->IsSlave( qn ) ) {
            queries_p_->SetDone4Search< true >( s->QNum() );
        }
    }
}

//------------------------------------------------------------------------------
void CBatch::MarkUPRes( CResult * s, CResult * e ) {
    for( ; s < e ; ++s ) {
        TQNum qn( s->QNum() ); 
        queries_p_->SetUPRes( qn );
        queries_p_->SetUPRes( queries_p_->GetMate( qn ) );
    }
}

//------------------------------------------------------------------------------
template< int search_mode >
void CBatch::RetainUnpairedResults( CResult * s, CResult * e )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    CResult * ss( s );

    for( ; s < e; ++s ) {
        TQNum qn( s->QNum() ), pqn( queries_p_->PrimaryQNum( qn ) );

        if( qn == pqn ) {
            if( queries_p_->DelResult< TScoring >( qn, *s ) ) {
                // CResult r( p_tmpres_mgr_.Save( s->GetRawLen() ) );
                CResult r( p_tmpres_mgr_->Save( s->GetRawLen() ) );
                r.Clone( *s );
            }
        }
    }

    for( s = ss; s < e; ++s ) {
        TQNum qn( s->QNum() ), pqn( queries_p_->PrimaryQNum( qn ) );
        if( qn == pqn ) queries_p_->AddResult< TScoring >( qn, *s );
    }
}

//------------------------------------------------------------------------------
void CBatch::ClearDone4Search( CResult * s, CResult * e )
{ for( ; s < e; ++s ) queries_p_->SetDone4Search< false >( s->QNum(), false ); }

//------------------------------------------------------------------------------
template< int search_mode >
bool CBatch::CombineUnpairedResults( 
        CResult * s, CResult * e, 
        size_t & n_found, size_t & n_left, size_t & n_cancelled )
{
    bool res( false );
    CResult * m( s );
    while( m != e && queries_p_->IsLeft( m->QNum() ) ) ++m;
    size_t r1( m - s ), r2( e - m );
    size_t ru( queries_p_->ResLimit() );
    RetainUnpairedResults< search_mode >( s, e );

    if( r1 == 0 || r2 == 0 ) {
        ++n_cancelled;
        UpdateCounts< search_mode >( s, e );
        ClearDone4Search( s, e );
        MarkDone( s, e );
    }
    else if( r1 < ru && r2 < ru ) {
        if( IdentifyPairs< search_mode >( s, r1, r2 ) == 0 ) { 
            ++n_left; 
            ClearDone4Search( s, e );
            res = true;
        }
        else {
            ++n_found;
            UpdateCounts< search_mode >( s, e );
            MarkDone( s, e );
            MarkUPRes( s, e );
        }
    }
    else { 
        ++n_left;
        ClearDone4Search( s, e );
        res = true;
    }

    return res;
}

//------------------------------------------------------------------------------
void CBatch::UpdateHistogramForSubject( 
        CResult * ls, CResult * le, CResult * rs, CResult * re, THistogram * h,
        bool * scv )
{
    for( CResult * i( ls ); i != le; ++i ) {
        for( CResult * j( rs ); j != re; ++j ) {
            if( CheckStrandConfig( i, j ) ) {
                int sc( GetStrandConfig( i, j ) );
                Uint4 d( (i->SOff( 0 ) < j->SOff( 0 )) ?
                            j->SOff( 0 ) - i->SOff( 0 ) :
                            i->SOff( 0 ) - j->SOff( 0 ) );

                if( d < ISD_MAX_INSERT && !scv[sc] ) {
                    scv[sc] = true;
                    ++h[sc][d];
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
void CBatch::UpdateHistogram( CResult * s, CResult * e, THistogram * h )
{
    CResult * m( s );
    while( m != e && queries_p_->IsLeft( m->QNum() ) ) ++m;
    if( m == s || m == e ) return;
    bool scv[4];
    std::fill( scv, scv + 4, false );

    for( CResult * i( s ), * j( m ); i != m; ) {
        TDBOrdId snum( i->SNum() );
        while( j != e && j->SNum() < snum ) ++j;
        if( j == e ) break;
        CResult * ie( i );
        while( ie != m && ie->SNum() == snum ) ++ie;

        if( j->SNum() == snum ) {
            CResult * je( j );
            while( je != e && je->SNum() == snum ) ++je;
            UpdateHistogramForSubject( i, ie, j, je, h, scv );
            j = je;
        }

        i = ie;
    }
}

//------------------------------------------------------------------------------
bool CBatch::ProcessSample( 
        TSeqSize & min, TSeqSize range, Uint8 total, const THistogram & h )
{
    Uint8 ctotal( 0 ), max_total( 0 );
    TSeqSize my_min, max;
    THistogram::const_iterator min_i( h.begin() ), max_i( min_i );

    while( max_i != h.end() ) {
        my_min = min_i->first;
        max = my_min + range;
        
        while( max_i != h.end() && max > max_i->first ) {
            ctotal += max_i->second;
            ++max_i;
        }

        if( ctotal >= ISD_ACCEPT_PER_SC_RATIO*total &&
                ctotal > max_total ) {
            min = my_min;
            max_total = ctotal;
        }

        ++min_i;
        max_i = min_i;
        ctotal = 0;
    }

    return max_total > 0;
}

//------------------------------------------------------------------------------
void CBatch::DumpHistogram( THistogram * h )
{
    M_TRACE( CTracer::INFO_LVL, "dumping histogram data" );
    std::ofstream os( init_data_.hist_fname.c_str() );

    if( !os ) {
        M_TRACE( CTracer::INFO_LVL, 
                "could not open " << init_data_.hist_fname );
        return;
    }

    for( int i( 0 ); i < 4; ++i ) {
        if( !h[i].empty() ) {
            os << "strand configuration " << i << std::endl;
                    
            for( THistogram::const_iterator j( h[i].begin() ); 
                    j != h[i].end(); ++j ) {
                os << j->first << ' ' << j->second << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------
bool CBatch::ProcessHistogram( THistogram * h )
{
    Uint8 totals[4], total( 0 );
    std::fill( totals, totals + 4, 0 );

    // Compute totals for each strand configuration, and overall total;
    //
    for( int i( 0 ); i < 4; ++i ) {
        for( THistogram::const_iterator j( h[i].begin() );
                j != h[i].end(); ++j ) {
            totals[i] += j->second;
        }

        total += totals[i];
    }

    // Select the one with ISD_ACCEPT_RATIO values.
    //
    int i( 0 );
    for( ; i < 4; ++i ) if( totals[i] >= ISD_ACCEPT_RATIO*total ) break;
    if( i == 4 ) return false;

    // Compute the candidate insert size and deviation
    //
    TSeqSize min( 0 ), range( 0 ), tmin( 0 ), 
             trange( ISD_MAX_RANGE ), f( trange );

    do {
        if( ProcessSample( tmin, trange, totals[i], h[i] ) ) {
            min = tmin;
            range = trange;
            f /= 2;
            trange -= f;
            if( trange < 128 ) break;
        }
        else if( range == 0 ) return false;
        else {
            f /= 2;
            trange += f;
        }
    }
    while( f > 1 );

    init_data_.discover_sep = false;
    range /= 2;
    init_data_.pair_distance = min + range;
    TSeqSize max_range( init_data_.pair_distance/10 + ISD_MIN_RANGE );
    max_range = std::max( max_range, (TSeqSize)init_data_.pair_fuzz );
    M_TRACE( CTracer::INFO_LVL,
             "candidate trace configuration: " << i << ";" );
    M_TRACE( CTracer::INFO_LVL,
             "candidate insert size: " << init_data_.pair_distance << ";" );
    M_TRACE( CTracer::INFO_LVL,
             "candidate insert size deviation: " << range );

    if( range > max_range ) {
        M_TRACE( CTracer::INFO_LVL,
                 "candidate deviation is over max of: " <<  max_range );
        return false;
    }

    init_data_.pair_fuzz = std::max( (Uint2)(range), init_data_.pair_fuzz );

    {
        init_data_.resconf_str = "0000";
        init_data_.resconf_str[i] = '1';
        init_data_.ipam_vec = ParseResConfStr( init_data_.resconf_str );
    }

    pass_init_data_.pair_distance = init_data_.pair_distance;
    pass_init_data_.pair_fuzz = init_data_.pair_fuzz;
    pass_init_data_.ipam_vec = init_data_.ipam_vec;

    M_TRACE( CTracer::INFO_LVL,
             "selected strand configuration: " << i );
    M_TRACE( CTracer::INFO_LVL,
             "selected insert size: " << init_data_.pair_distance );
    M_TRACE( CTracer::INFO_LVL,
             "selected insert range: " << init_data_.pair_fuzz );
    return true;
}

//------------------------------------------------------------------------------
template< int search_mode >
bool CBatch::DiscoverInsertSize( void )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    M_TRACE( CTracer::INFO_LVL, "attempting to discover insert size" );
    queries_p_->ReadBackQueryData( tmp_store_ );
    queries_p_->CleanUpQueryData();
    queries_p_->ClearMarks();
    queries_p_->FreeQueryData();
    queries_p_->StartPostProcess();
    queries_p_->DistributeUnpairedCounts< TScoring, false >();
    CMemoryManager * mem_mgr_p( pass_init_data_.mem_mgr_p );

    {
        THistogram histogram[4];

        SResultsHolder rh( *mem_mgr_p );
        char * res_data_start( (char *)rh.Get() ), 
             * res_data_end( res_data_start );
        CResult * res_end( (CResult *)rh.Get() + rh.Size()/sizeof( CResult ) ),
                * res_start( res_end );
        size_t max_results( 0 );
        size_t res_limit( ISD_RES_LIMIT );

        {
            size_t res_len( CResult::EstimateLen( 1, pass_init_data_.n_err ) );
            size_t free_space( (char *)res_end - res_data_start );
            max_results = free_space/(res_len + sizeof( CResult ));
        }

        std::pair< TQNum, TQNum > bounds = 
            std::make_pair< TQNum, TQNum >( 0UL, 0UL );
    
        while( bounds.second < queries_p_->size() ) {
            bounds = ComputeQNumBounds< search_mode, true >( 
                    bounds.second, max_results );
            M_TRACE( CTracer::INFO_LVL, 
                     "processing results for queries " << bounds.first <<
                     " -- " << bounds.second );
            M_TRACE( CTracer::INFO_LVL, "loading results" );
            // u_tmpres_mgr_.LoadInit();
            u_tmpres_mgr_->LoadInit();
            size_t n_res = 0;

            while( n_res < max_results ) {
                // *(--res_start) = u_tmpres_mgr_.Load();
                *(--res_start) = u_tmpres_mgr_->Load();
                if( res_start->Empty() ) { ++res_start; break; }
                TQNum qn( res_start->QNum() );

                if( queries_p_->NRes< TScoring, false >( qn ) >= res_limit ) {
                    ++res_start; continue;
                }

                if( queries_p_->IsUnique( qn ) ) {
                    if( queries_p_->HasHigherRank< TScoring >( 
                                qn, *res_start ) ) {
                        ++res_start; continue;
                    }

                    res_start->Copy( res_data_end );
                    res_data_end += res_start->GetRawLen();
                    ++n_res;
                }
                else {
                    TQNum * dis( queries_p_->DupStart( qn ) ),
                          * die( queries_p_->DupEnd( qn ) );

                    for( TQNum * di = dis; di != die; ++di ) {
                        if( *di >= bounds.first && *di < bounds.second ) {
                            if( queries_p_->HasHigherRank< TScoring >( 
                                        *di, *res_start ) ) {
                                break;
                            }
                        }

                        *(res_start - 1) = *res_start;
                        res_start->Copy( res_data_end );
                        res_data_end += res_start->GetRawLen();
                        res_start->SetQNum( *di );
                        ++n_res;
                        --res_start;
                    }

                    ++res_start;
                }
            }

            M_TRACE( CTracer::INFO_LVL, "loaded " << n_res << " results" );
            std::stable_sort( 
                    res_start, res_end, CResult::SHLCompare( &seqstore_ ) );
            M_TRACE( CTracer::INFO_LVL, "results sorted" );
            CResult * s( res_start ), * e( s ), * rend( s + n_res );

            while( s != rend ) {
                while( e != rend && 
                        e->QOrdId( start_qid_, true ) ==
                            s->QOrdId( start_qid_, true ) ) {
                    ++e;
                }

                UpdateHistogram( s, e, histogram );
                s = e;
            }

            res_data_end = res_data_start;
            res_start = res_end;
        }

        DumpHistogram( histogram );

        if( !ProcessHistogram( histogram ) ) {
            M_TRACE( CTracer::INFO_LVL, "insert size discovery failed" );
            return false;
        }

        return true;
    }
}

//------------------------------------------------------------------------------
template< int search_mode >
bool CBatch::InterProcess( void )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    size_t n_found( 0 ), n_left( 0 ), n_cancelled( 0 );

    if( init_data_.discover_sep ) {
        bool cont( DiscoverInsertSize< search_mode >() );
        cont = cont && !init_data_.discover_sep_stop;

        if( !cont ) {
            queries_p_->SetResLimit< TScoring >( res_limit_ );
            return false;
        }

        M_TRACE( CTracer::INFO_LVL, "analyzing single results" );
    }
    else {
        M_TRACE( CTracer::INFO_LVL, "analyzing single results" );
        queries_p_->ReadBackQueryData( tmp_store_ );
        queries_p_->CleanUpQueryData();
        queries_p_->ClearMarks();
        queries_p_->FreeQueryData();
        queries_p_->StartPostProcess();
        queries_p_->DistributeUnpairedCounts< TScoring, false >();
    }

    n_cancelled += queries_p_->MarkDoneNotFound< TScoring >();
    CMemoryManager * mem_mgr_p( pass_init_data_.mem_mgr_p );
    std::ofstream os;

    if( !paired_log_.empty() ) {
        try {
            os.open( paired_log_.c_str() );
        }
        catch( ... ) {
            M_TRACE( CTracer::WARNING_LVL, 
                     "can not open paired log file " << paired_log_ <<
                     " for writing" );
            paired_log_ = "";
        }
    }

    {
        SResultsHolder rh( *mem_mgr_p );
        char * res_data_start( (char *)rh.Get() ), 
             * res_data_end( res_data_start );
        CResult * res_end( (CResult *)rh.Get() + rh.Size()/sizeof( CResult ) ),
                * res_start( res_end );
        size_t max_results( 0 );

        {
            size_t res_len( CResult::EstimateLen( 1, pass_init_data_.n_err ) );
            size_t free_space( (char *)res_end - res_data_start );
            max_results = free_space/(res_len + sizeof( CResult ));
        }

        std::pair< TQNum, TQNum > bounds = 
            std::make_pair< TQNum, TQNum >( 0UL, 0UL );
    
        while( bounds.second < queries_p_->size() ) {
            bounds = ComputeQNumBounds< search_mode, true >( 
                    bounds.second, max_results );
            M_TRACE( CTracer::INFO_LVL, 
                     "processing results for queries " << bounds.first <<
                     " -- " << bounds.second );
            M_TRACE( CTracer::INFO_LVL, "loading results" );
            // u_tmpres_mgr_.LoadInit();
            u_tmpres_mgr_->LoadInit();
            size_t n_res = 0;

            while( n_res < max_results ) {
                // *(--res_start) = u_tmpres_mgr_.Load();
                *(--res_start) = u_tmpres_mgr_->Load();
                if( res_start->Empty() ) { ++res_start; break; }
                TQNum qn( res_start->QNum() );

                if( queries_p_->IsUnique( qn ) ) {
                    if( queries_p_->HasHigherRank< TScoring >( 
                                qn, *res_start ) ) {
                        ++res_start; continue;
                    }

                    if( qn >= bounds.first && qn < bounds.second &&
                            queries_p_->DelResult< TScoring >( 
                                qn, *res_start ) ) {
                        res_start->Copy( res_data_end );
                        res_data_end += res_start->GetRawLen();
                        ++n_res;
                    }
                    else ++res_start;
                }
                else {
                    TQNum * dis( queries_p_->DupStart( qn ) ),
                          * die( queries_p_->DupEnd( qn ) );

                    for( TQNum * di = dis; di != die; ++di ) {
                        if( *di >= bounds.first && *di < bounds.second ) {
                            if( !queries_p_->HasPairedResult< TScoring >( 
                                        *di ) ) {
                                if( queries_p_->HasHigherRank< TScoring >( 
                                            *di, *res_start ) ) {
                                    break;
                                }

                                if( queries_p_->DelResult< TScoring >( 
                                            *di, *res_start ) ) {
                                    *(res_start - 1) = *res_start;
                                    res_start->Copy( res_data_end );
                                    res_data_end += res_start->GetRawLen();
                                    res_start->SetQNum( *di );
                                    ++n_res;
                                    --res_start;
                                }
                            }
                        }
                    }

                    ++res_start;
                }
            }

            // u_tmpres_mgr_.LoadFinal();
            u_tmpres_mgr_->LoadFinal();
            M_TRACE( CTracer::INFO_LVL, "loaded " << n_res << " results" );
            std::stable_sort( 
                    res_start, res_end, CResult::SHLCompare( &seqstore_ ) );
            M_TRACE( CTracer::INFO_LVL, "results sorted" );
            CResult * s( res_start ), * e( s ), * rend( s + n_res );

            while( s != rend ) {
                while( e != rend && 
                        e->QOrdId( start_qid_, true ) ==
                            s->QOrdId( start_qid_, true ) ) {
                    queries_p_->AddResult< TScoring >( e->QNum(), *e );
                    ++e;
                }

                bool left( CombineUnpairedResults< search_mode >( 
                            s, e, n_found, n_left, n_cancelled ) );

                if( !paired_log_.empty() && left ) {
                    os << s->QOrdId( start_qid_, true ) << '\n';
                }

                s = e;
            }

            res_data_end = res_data_start;
            res_start = res_end;
        }
    }

    if( !paired_log_.empty() ) os.flush();

    queries_p_->MarkDoneSlaves();
    queries_p_->SetResLimit< TScoring >( res_limit_ );
    queries_p_->ReadPrimaryQueryData( tmp_store_ );
    queries_p_->ClearHashUseInfo();
    queries_p_->ClearMarks();

    M_TRACE( CTracer::INFO_LVL, 
             "number of queries with paired results: " << n_found );
    M_TRACE( CTracer::INFO_LVL,
             "number of queries left for paired search: " << n_left );
    M_TRACE( CTracer::INFO_LVL,
             "number of cancelled queries: " << n_cancelled );

    return true;
}

END_NS( srprism )
END_STD_SCOPES

