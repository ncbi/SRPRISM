/*  $Id: batch_priv.hpp 539744 2017-06-27 13:06:13Z morgulis $
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

#include <utility>
#include <memory>
#include <algorithm>
#include <cstdlib>

#ifndef NCBI_CPP_TK

#include "sidmap.hpp"
#include "query_store.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/srprism/sidmap.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_store.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
namespace {
    template< bool paired > struct SQueryGroupSize {};
    template<> struct SQueryGroupSize< false > { static const int VALUE = 1; };
    template<> struct SQueryGroupSize< true >  { static const int VALUE = 2; };
}

//------------------------------------------------------------------------------
template< int hash, bool paired > inline void CBatch::RunPass( 
        bool skip_good, bool skip_bad )
{
    if( search_mode_ == SSearchMode::PARTIAL ) {
        try {
            CSearchPass< SSearchMode::PARTIAL, hash, paired > pass( 
                    pass_init_data_ );
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            pass.Run();
        }
        catch( typename CSearchPass< 
                SSearchMode::PARTIAL, hash, paired >::CException & e ) {
            if( e.ErrorCode() == CSearchPass< 
                        SSearchMode::PARTIAL, hash, paired 
                    >::CException::PASS_SKIP ) {
                M_TRACE( CTracer::INFO_LVL, e.what() );
            }
            else throw;
        }
    }
    else if( search_mode_ == SSearchMode::DEFAULT ) {
        try {
            CSearchPass< SSearchMode::DEFAULT, hash, paired > pass( 
                    pass_init_data_ );
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            pass.Run();
        }
        catch( typename CSearchPass< 
                SSearchMode::DEFAULT, hash, paired >::CException & e ) {
            if( e.ErrorCode() == CSearchPass< 
                    SSearchMode::DEFAULT, 
                    hash, 
                    paired >::CException::PASS_SKIP ) {
                M_TRACE( CTracer::INFO_LVL, e.what() );
            }
            else throw;
        }
    }
    else if( search_mode_ == SSearchMode::SUM_ERR ) {
        try {
            CSearchPass< SSearchMode::SUM_ERR, hash, paired > pass( 
                    pass_init_data_ );
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            pass.Run();
        }
        catch( typename CSearchPass< 
                SSearchMode::SUM_ERR, hash, paired >::CException & e ) {
            if( e.ErrorCode() == CSearchPass< 
                    SSearchMode::SUM_ERR, 
                    hash, 
                    paired >::CException::PASS_SKIP ) {
                M_TRACE( CTracer::INFO_LVL, e.what() );
            }
            else throw;
        }
    }
    else if( search_mode_ == SSearchMode::BOUND_ERR ) {
        try {
            CSearchPass< SSearchMode::BOUND_ERR, hash, paired > pass( 
                    pass_init_data_ );
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            if( skip_good ) pass.SkipGood();
            if( skip_bad ) pass.SkipBad();
            pass.Run();
        }
        catch( typename CSearchPass< 
                SSearchMode::BOUND_ERR, hash, paired >::CException & e ) {
            if( e.ErrorCode() == CSearchPass< 
                    SSearchMode::BOUND_ERR, 
                    hash, 
                    paired >::CException::PASS_SKIP ) {
                M_TRACE( CTracer::INFO_LVL, e.what() );
            }
            else throw;
        }
    }
    else SRPRISM_ASSERT( false );
}

//------------------------------------------------------------------------------
namespace {
    struct SResultsHolder
    {
        SResultsHolder( CMemoryManager & mem_mgr )
            : mem_mgr_( mem_mgr ), ptr_( 0 )
        {
            sz_ = mem_mgr_.GetFreeSpace();
            ptr_ = mem_mgr_.Allocate( mem_mgr_.GetFreeSpace() );
        }

        ~SResultsHolder(void) { mem_mgr_.Free( ptr_ ); }

        void * Get(void) const { return ptr_; }
        size_t Size( void ) const { return sz_; }

        private:

            CMemoryManager & mem_mgr_;
            void * ptr_;
            size_t sz_;
    };

    Sint2 Quality( size_t n_res, size_t max_res )
    { return n_res < max_res ? std::max( 100/n_res, 1UL ) : 0; }

    struct SALCounts
    {
        SALCounts( const CSeqStore & ss )
            : counts( ss.NSeq(), std::make_pair( (TQueryOrdId)0, (size_t)0 ) ),
              rcounts( ss.NSeq(), std::make_pair( (TQueryOrdId)0, (size_t)0 ) )
        {}

        void Add( TDBOrdId sid, TQueryOrdId qid )
        {
            if( qid != counts[sid].first ) {
                counts[sid] = std::make_pair( qid, (size_t)0 );
            }

            ++counts[sid].second;
        }

        void AddR( TDBOrdId sid, TQueryOrdId qid )
        {
            if( qid != rcounts[sid].first ) {
                rcounts[sid] = std::make_pair( qid, (size_t)0 );
            }

            ++rcounts[sid].second;
        }

        size_t GetCount( TDBOrdId sid, TQueryOrdId qid ) const
        { return (qid == counts[sid].first) ? counts[sid].second : 0; }

        size_t GetRCount( TDBOrdId sid, TQueryOrdId qid ) const
        { return (qid == rcounts[sid].first) ? rcounts[sid].second : 0; }


        private:

            typedef std::vector< std::pair< TQueryOrdId, size_t > > TData;
            TData counts, rcounts;
    };
}

//------------------------------------------------------------------------------
namespace {

template< typename t_scoring_sys >
struct SCompareForOutput {
    typedef t_scoring_sys TScoring;

    SCompareForOutput( CSeqStore * s ) : hlc( s ) {}

    bool operator()( CResult const * l, CResult const * r ) const { 
        CResult const & ll( *l ), & rr( *r );

        if( !hlc( ll, rr ) ) {
            if( !hlc( rr, ll ) ) {
                if( !CResult::SLLCompare()( ll, rr ) ) {
                    if( !CResult::SLLCompare()( rr, ll ) ) {
                        return CResult::CompareLevels< TScoring >()( ll, rr );
                    }
                    else return false;
                }
                else return true;
            }
            else return false;
        }
        else return true;
    }

private:

    CResult::SHLCompare hlc;
};

}

//------------------------------------------------------------------------------
template< int search_mode, bool paired > 
CResult * CBatch::RemoveDuplicatesForSubject( 
        CResult * s, CResult * e, CResult * o ) {
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    if( s == e ) return o;
    common::Uint1 nerr( pass_init_data_.n_err );
    CResult * h( s ), * b( h + 1 ), * w( h + 1 ), * t( h );

    while( h != e ) {
        while( w != e && 
               abs( (Sint8)(h->SOff( 0 )) - w->SOff( 0 ) ) <= nerr ) {
            if( !paired || !s->Paired() || 
                    abs( (Sint8)(h->SOff( 1 )) - w->SOff( 1 ) ) <= nerr ) {
                std::swap( *b, *w );
                if( CResult::CompareLevels< TScoring >()( *b, *h ) ) t = b;
                ++b;
            }

            ++w;
        }

        *o = *t;
        if( t != h ) { h = t; w = b; } else { h = t = b++; w = b; ++o; }
    }

    return o;
}

//------------------------------------------------------------------------------
template< int search_mode, bool paired > 
CResult * CBatch::RemoveDuplicates( CResult * s, CResult * e ) {
    CResult * o( s );

    while( s != e ) {
        CResult * se( s + 1 );
        for( ; se != e && s->SNum() == se->SNum(); ++se );
        std::stable_sort( s, se, CResult::SLLCompare() );

        while( s != se ) {
            CResult * sse( s + 1 );
            for( ; sse != e && CResult::SLLCompare::IsEquiv( *s, *sse ); 
                    ++sse );
            o = RemoveDuplicatesForSubject< search_mode, paired >( s, sse, o );
            s = sse;
        }
    }

    return o;
}

//------------------------------------------------------------------------------
template< int search_mode, bool paired >
void CBatch::PostProcess(void)
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    M_TRACE( CTracer::INFO_LVL, "post processing results" );
    queries_p_->ReadBackQueryData( tmp_store_ );
    queries_p_->CleanUpQueryData();
    queries_p_->ClearMarks();
    queries_p_->FreeQueryData();
    queries_p_->StartPostProcess();
    queries_p_->DistributeUnpairedCounts< TScoring, paired >();
    if( !paired ) queries_p_->ClearHashUseInfo();
    CMemoryManager * mem_mgr_p( pass_init_data_.mem_mgr_p );
    SResultsHolder rh( *mem_mgr_p );
    char * res_data_start( (char *)rh.Get() ), 
         * res_data_end( res_data_start );
    CResult * res_end( (CResult *)rh.Get() + rh.Size()/sizeof( CResult ) ),
            * res_start( res_end );
    size_t max_results( 0 );

    {
        size_t res_len( 
                !paired ? CResult::EstimateLen( 1, pass_init_data_.n_err )
                        : CResult::EstimateLen( 
                           2, pass_init_data_.n_err, pass_init_data_.n_err ) );
        size_t free_space( (char *)res_end - res_data_start );
        max_results = free_space/(res_len + sizeof( CResult ));
        M_TRACE( CTracer::INFO_LVL, "space for " << max_results << " results" );
    }

    std::pair< TQNum, TQNum > bounds = 
        std::make_pair< TQNum, TQNum >( 0UL, 0UL );
    SALCounts al_counts( seqstore_ );
    typedef std::vector< TDBOrdId > TALIds;
    typedef TALIds::const_iterator TALIdsIter;
    TALIds al_ids;

    out_p_->SetUpQueryInfo( 
            queries_p_.get(), paired ? start_qid_/2 : start_qid_ );

    while( bounds.second < queries_p_->size() ) {
        bounds = ComputeQNumBounds< search_mode, paired >( 
                bounds.second, max_results );
        M_TRACE( CTracer::INFO_LVL, 
                 "processing results for queries " << bounds.first <<
                 " -- " << bounds.second );
        tmpres_mgr_p_->LoadInit();
        size_t n_res = 0;

        while( n_res < max_results ) {
            *(--res_start) = tmpres_mgr_p_->Load();
            if( res_start->Empty() ) { ++res_start; break; }

            if( paired && !res_start->Paired() ) {
                if( !seqstore_.CheckRegion(
                            seqstore_.EncodePos( std::make_pair(
                                    res_start->SNum(),
                                    res_start->SOff( 0 ) ) ),
                                    // res_start->SOff( 0 ) + res_start->GetLeftOffset( 0 ) ) ),
                            res_start->GetAlignLen( 0 ) - 
                                res_start->GetNIns( 0 ) ) ) {
                    ++res_start;
                    continue;
                }
            }

            TQNum qn( res_start->QNum() );

            if( res_start->Paired() || queries_p_->IsUnique( qn ) ) {
                if( queries_p_->HasHigherRank< TScoring >( qn, *res_start ) ) {
                    ++res_start; continue;
                }

                if( qn >= bounds.first && qn < bounds.second &&
                        queries_p_->DelResult< TScoring >( qn, *res_start ) ) {
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
                        if( !queries_p_->HasPairedResult< TScoring >( *di ) ) {
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

        tmpres_mgr_p_->LoadFinal();
        M_TRACE( CTracer::INFO_LVL, "loaded " << n_res << " results" );
        std::stable_sort( res_start, res_end, CResult::SHLCompare( &seqstore_ ) );
        M_TRACE( CTracer::INFO_LVL, "results sorted" );

        if( queries_p_->QueriesReversed() ) {
            for( CResult * i( res_start ); i != res_end; ++i ) {
                i->ReverseStrands();
            }

            M_TRACE( CTracer::INFO_LVL, "strands reversed" );
        }

        typedef std::vector< CResult * > TResPtrSet;
        TResPtrSet results;

        for( CResult * rstart( res_start ), 
                     * rend( rstart + n_res ); rstart != rend; ) {
            CResult * qrend( rstart );
            
            {
                TQueryOrdId qid( qrend->QOrdId( start_qid_, paired ) );
                do { ++qrend; } 
                while( qrend != rend && 
                       qid == qrend->QOrdId( start_qid_, paired ) );
            }

            results.clear();

            while( rstart != qrend ) {
                size_t n_ref( 0 );
                TQNum qn( rstart->QNum() );
                CResult * rrend( rstart );
                do { ++rrend; } while( rrend != qrend && rrend->QNum() == qn );
                size_t n_res( rrend - rstart );

                // zero quality: number of results found is greater than
                // internal limit; in this case quality of all results
                // is set to 0 and no duplication removal is performed
                //
                bool zero_quality( n_res >= res_limit_ );

                CResult * rdrend( 
                        zero_quality ? rrend 
                                     : RemoveDuplicates< search_mode, paired >( 
                                            rstart, rrend ) );

                if( !zero_quality ) {
                    // count the number of results on primary that overlap 
                    // alternate loci regions
                    //
                    for( CResult * r( rstart ); r != rdrend; ++r ) {
                        TDBOrdId sid( r->SNum() ), 
                                 rsid( seqstore_.GetRefOId( sid ) );
                        al_counts.Add( sid, qn );

                        if( sid == rsid ) { // primary seq
                            ++n_ref;
                            al_ids.clear();
                            seqstore_.GetInsideList( 
                                        rsid, r->SOff( 0 ), al_ids );
                            seqstore_.GetInsideList( 
                                    rsid, 
                                    r->SOff( 0 ) + r->GetAlignLen( 0 ) - 
                                        r->GetNIns( 0 ), 
                                    al_ids );

                            if( r->Paired() ) {
                                seqstore_.GetInsideList( 
                                        rsid, r->SOff( 1 ), al_ids );
                                seqstore_.GetInsideList( 
                                        rsid, r->SOff( 1 ) + 
                                            r->GetAlignLen( 1 ) - 
                                            r->GetNIns( 1 ), 
                                        al_ids );
                            }

                            std::stable_sort( al_ids.begin(), al_ids.end() );
                            TALIdsIter e( std::unique( 
                                        al_ids.begin(), al_ids.end() ) );

                            for( TALIdsIter ali( al_ids.begin() ); 
                                    ali != e; ++ali ) {
                                al_counts.AddR( *ali, qn );
                            }
                        }
                    }
                }

                size_t final_n_res( 
                        std::min( (Uint4)n_res, final_res_limit_ - 1 ) );
                CResult * sstart( rstart ), * send( sstart );

                while( sstart != rdrend ) {
                    // select a group of results with the same subject id;
                    // stop when primary sequence is reached
                    //
                    TDBOrdId sid( sstart->SNum() );
                    if( sid == seqstore_.GetRefOId( sid ) ) break; // primary
                    for( ; send != rdrend && send->SNum() == sid; ++send );
                    Sint2 qual( 0 );

                    if( !zero_quality ) {
                        size_t nr( n_ref );
                        nr += al_counts.GetCount( sid, qn );
                        nr -= al_counts.GetRCount( sid, qn );
                        nr = std::min( nr, (size_t)final_res_limit_ );
                        qual = Quality( nr, final_res_limit_ );
                    }

                    std::stable_sort( sstart, send, 
                               CResult::CompareLevels< TScoring >() );

                    for( size_t i( 0 ); i < final_n_res && sstart != send; 
                            ++i, ++sstart ) {
                        sstart->SetQuality( qual );
                        results.push_back( sstart );
                    }

                    sstart = send;
                }

                // processing results for primary sequences
                //
                if( !zero_quality && n_ref <= final_n_res ) {
                    for( size_t i( 0 ); i < final_n_res && sstart != rdrend; 
                            ++sstart, ++i ) {
                        sstart->SetQuality( 
                                Quality( n_ref, final_res_limit_ ) );
                        results.push_back( sstart );
                    }
                }
                else {
                    // we have more results on primary than needed,
                    // so diversify subjects at the worst level
                    //
                    std::stable_sort( sstart, rdrend, 
                               CResult::CompareLevels< TScoring >() );
                    CResult * lstart( sstart ), * lend( lstart );

                    while( lstart != lend && 
                                (size_t)(lend - sstart) < final_n_res ) {
                        /*
                            Find the range of results (lstart,lend) with the 
                            same level as lstart. 
                            If lend - sstart > final_n_res, then only keep 
                            final_n_res - (lend - sstart) from (lstart, lend)
                            with as many different subjects as possible;
                            otherwise keep all of (lstart, lend).
                        */
                        while( lend != rdrend &&
                                TScoring::HaveEqualLevels( *lstart, *lend ) ) {
                            ++lend;
                        }

                        if( (size_t)(lend - sstart) > final_n_res ) {
                            DiversifySubjects( lstart, lend );
                            break;
                        }

                        lstart = lend;
                    }

                    lend = sstart + std::min( 
                            final_n_res, (size_t)(rdrend - sstart) );

                    for( ; sstart != lend; ++sstart ) {
                        sstart->SetQuality( 0 );
                        results.push_back( sstart );
                    }
                }

                for( ; rstart != rdrend; ++rstart ) {
                    queries_p_->AddResult< TScoring >( qn, *rstart );
                }

                queries_p_->SetMark( qn, false );
                rstart = rrend;
            }

            //
            // sort for output
            //
            std::stable_sort( results.begin(), results.end(), 
                              SCompareForOutput< TScoring >( &seqstore_ ) );

            // compute if the results are guaranteed best
            //
            int pg[2] = { 0, 0 };

            if( search_mode_ == SSearchMode::DEFAULT ) {
                if( !results.empty() ) {
                    CResult * r( *results.begin() );

                    if( paired ) {
                        TQNum qn1( r->QNum() ), qn2;

                        if( !queries_p_->IsLeft( qn1 ) ) {
                            qn1 = queries_p_->GetMate( qn1 );
                        }

                        qn2 = queries_p_->GetMate( qn1 );
                        TQNum qm( qn1 ), qs( qn2 );
                        if( queries_p_->IsSlave( qm ) ) std::swap( qm, qs );

                        if( r->Paired() ) {
                            if( queries_p_->IsUPRes( qn1 ) ) {
                                if( ProvidesGuarantee( qn1, r->NErr( 0 ) ) &&
                                    ProvidesGuarantee( qn2, r->NErr( 1 ) ) )
                                { 
                                    pg[0] = pg[1] = 1;
                                }
                                else pg[0] = pg[1] = 0;
                            }
                            else {
                                int nerr( std::max( 
                                                r->NErr( 0 ), r->NErr( 1 ) ) );
                                pg[0] = pg[1] = ProvidesGuarantee( qm, nerr );
                            }
                        }
                        else {
                            TResPtrSet::const_iterator i( results.begin() ),
                                                       j( i );

                            while( j != results.end() && (*j)->QNum() == qn1 ) {
                                ++j;
                            }

                            if( i != j && j != results.end() ) {
                                if( !queries_p_->HasRepHashes( qm ) ) {
                                    if( qm == qn1 ) {
                                        pg[0] = 1;
                                        pg[1] = ProvidesGuarantee( 
                                                qs, (*j)->NErr( 0 ) );
                                    }
                                    else {
                                        pg[1] = 1;
                                        pg[0] = ProvidesGuarantee( 
                                                qs, (*i)->NErr( 0 ) );
                                    }
                                }
                            }
                            else if( i != j ) {
                                if( !queries_p_->HasRepHashes( qn2 ) ) {
                                    pg[1] = 1;
                                    pg[0] = ProvidesGuarantee( 
                                            qn1, (*i)->NErr( 0 ) );
                                }
                            }
                            else if( j != results.end() ) {
                                if( !queries_p_->HasRepHashes( qn1 ) ) {
                                    pg[0] = 1;
                                    pg[1] = ProvidesGuarantee( 
                                            qn2, (*j)->NErr( 0 ) );
                                }
                            }
                        }
                    }
                    else {
                        pg[0] = ProvidesGuarantee( 
                                r->QNum(), r->NErr( 0 ) );
                    }
                }
            }

            out_p_->ResultsOut( results, start_qid_, pg );
        }

        res_data_end = res_data_start;
        res_start = res_end;
    }

    out_p_->FinalizeBatch();
}

//------------------------------------------------------------------------------
template< int search_mode, bool paired >
std::pair< TQNum, TQNum > CBatch::ComputeQNumBounds( 
        TQNum start, size_t max_res ) const
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    static const int GROUP_SIZE = SQueryGroupSize< paired >::VALUE;
    size_t total = 0;
    TQNum end = start;
    bool have_res( false );

    while( end < queries_p_->size() ) {
        size_t group_res = 0;
        int i = 0;

        for( ; end + i < queries_p_->size() && i < GROUP_SIZE; ++i ) {
            group_res += queries_p_->NRes< TScoring, paired >( end + i );
        }

        if( total + group_res > max_res ) { have_res = true; break; }
        end += i;
        total += group_res;
    }

    if( total == 0 && have_res ) {
        M_THROW( CException, POSTPROCESS, "out of memory" );
    }

    return std::make_pair( start, end );
}

END_NS( srprism )
END_STD_SCOPES

