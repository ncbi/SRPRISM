/*  $Id: search_pass_priv.hpp 639115 2021-10-13 15:24:22Z morgulis $
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
 * File Description: one pass of a batch task
 *
 */

#include <iostream>

#include <cassert>
#include <functional>
#include <algorithm>

#ifndef NCBI_CPP_TK

#include <common/trace.hpp>
#include <common/binfile.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/result.hpp>

#include <srprism/query_data.hpp>
#include <srprism/query_store.hpp>
#include <srprism/inplace_align.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/trace.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/binfile.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>

#include <../src/internal/align_toolbox/srprism/lib/srprism/query_data.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_store.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqiter.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/inplace_align.hpp>

#include <../src/internal/align_toolbox/srprism/lib/srprism/stat.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
template< int search_mode > 
const char * CSearchPass_Base< search_mode >::PQDUMP_NAME = "pqdump";

//------------------------------------------------------------------------------
template< int search_mode >
CSearchPass_Base< search_mode >::CSearchPass_Base( 
        const CSearchPassDef::SInitData & init_data )
    : mem_mgr_( *init_data.mem_mgr_p ),
      tmp_res_mgr_( *init_data.tmpres_mgr_p ),
      seqstore_( *init_data.seqstore_p ),
      queries_( *init_data.queries_p ),
      tmp_store_( *init_data.tmp_store_p ),
      main_ma_( init_data.n_err, queries_.MaxQueryLen(), 1 ),
      ip_ma_( init_data.n_err, queries_.MaxQueryLen(), 1 ),
      idx_( nullptr ),
      idx_basename_( init_data.index_basename ),
      res_limit_( init_data.res_limit ),
      repeat_threshold_( init_data.repeat_threshold ),
      pair_distance_( init_data.pair_distance ), 
      pair_fuzz_( init_data.pair_fuzz ),
      n_err_( init_data.n_err ),
      end_pass_( false ),
      paired_search_( init_data.paired_search ),
      pass_stats_( init_data.search_stats ),
      bnf_( pass_stats_.n_filter, init_data.randomize ),
      randomize_( init_data.randomize ),
      random_seed_( init_data.random_seed )
{
    // if( random_seed_ ) srandom( time( 0 ) );
    if( random_seed_ ) srand( time( 0 ) );
    else srand( 1 );
    // else srandom( 1 );

    ext_data_table_.reserve( QEXT_TABLE_SIZE );
}

//------------------------------------------------------------------------------
namespace {
    // get max target number of errors depending on the query type
    //
    inline int MaxQueryErrors( bool is_short )
    { return is_short ? MAX_SHORT_ERR : MAX_ERR; }
}

//------------------------------------------------------------------------------
template< int search_mode, bool paired >
void CSearchPass_ByHash< search_mode, HASH_BLOWUP, paired >::NextSubPass( void )
{
    this->queries_.StartUpdate();
    start_idx_ = end_idx_;
    this->queries_.CleanUpQueryData();
    NextSubPass_Ext();
    this->queries_.SortQueries();
    this->queries_.SetupCrossLinks();
}

//------------------------------------------------------------------------------
namespace {
    inline int MinErrForBlowUp( const CQueryData & q )
    { return (q.IsShort() ? 1 : 2); }
}

//------------------------------------------------------------------------------
template< int search_mode, bool paired > void 
CSearchPass_ByHash< search_mode, HASH_BLOWUP, paired >::NextSubPass_Ext( void )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;

    if( paired ) {
        this->queries_.ClearMarks();

        CQueryData q( 0, 0, 0 );    // query data read from temp file
        size_t curr_idx( 0 ),       // index of the query currently being 
                                    //      processed
               expanded( 0 ),       // total number of expanded queries in this
                                    //      subpass
               extra( 0 );          // extra space needed for mates of expanded 
                                    //      queries

        std::vector< TQNum > marked_list;

        {
            CReadBinFile in( 
                    this->tmp_store_.Register( CQueryStore::QDUMP_NAME ) );

            while( in.Read( (char *)&q, sizeof( CQueryData ), true ) ) {
                size_t n_words( 3 + ((q.Len() - 1)>>WSHIFT) );
                n_words <<= 1;
                size_t extra_local( 0 );
                marked_list.clear();

                // reading raw query data
                //
                if( in.Read( 
                            (char *)&wbuf_[0], 
                            n_words*sizeof( TWord ), true ) == 0 ) {
                    SRPRISM_ASSERT( false );
                }

                // expand query if it belongs to the current subpass
                //
                if( curr_idx >= start_idx_ ) {
                    size_t n_q( 0 ); // number of query data entries resulting
                                     //     from expansion of one query
                    TQNum qn( q.QNum() );
                    q.SetRawData( &wbuf_[0] ); // link the query data to the raw 
                                               // sequence data

                    // if query does not have to be expanded, we still load it
                    // but mark as ignored, so it does not participate in the
                    // search
                    //
                    if( q.IsSALong() ||
                            TBase::template IgnoreForPass< true >( 
                                q, &this->queries_ ) ||
                            this->queries_.template GroupMaxErr< 
                                TScoring, true >( q.QNum() ) 
                                    < MinErrForBlowUp( q ) ||
                            this->n_err_ < MinErrForBlowUp( q ) ) {
                        q.SetIgnored( true );

                        if( !this->queries_.AddQueryData( 
                                &wbuf_[0], n_words, 
                                &q, 1, extra + extra_local ) ) {
                            break;
                        }
    
                        ++curr_idx; continue;
                    }

                    qcache_.clear();
                    BlowUp( q, n_q );

                    // add mates of 'q' and of its duplicates to the list of
                    // queries to load and adjust the needed extra space
                    //
                    if( this->queries_.IsUnique( qn ) ) {
                        TQNum pqn( this->queries_.PrimaryQNum( 
                                    this->queries_.GetMate( qn ) ) );
                        SetupMateMark( pqn, extra_local, marked_list );
                    }
                    else {
                        TQNum * ds( this->queries_.DupStart( qn ) ),
                              * de( this->queries_.DupEnd( qn ) );

                        for( ; ds != de; ++ds ) {
                            TQNum pqn( this->queries_.PrimaryQNum(
                                        this->queries_.GetMate( *ds ) ) );
                            SetupMateMark( 
                                    pqn, extra_local, marked_list );
                        }
                    }

                    if( n_q == 0 ) { // can happen for ambiguous queries
                        q.SetIgnored( true );

                        if( !this->queries_.AddQueryData(
                                    &wbuf_[0], n_words, &q, 1, 
                                    extra + extra_local ) ) {
                            break;
                        }
                        else {
                            for( std::vector< TQNum >::const_iterator i( 
                                        marked_list.begin() ); 
                                    i != marked_list.end(); ++i ) {
                                this->queries_.SetMark( *i, true );
                                extra += extra_local;
                            }

                            ++expanded;
                            ++curr_idx;
                            continue;
                        }
                    }

                    // try to append expanded query to the query store;
                    // if successful, mark all queries in the marked list
                    //
                    if( !this->queries_.AddQueryData( 
                                &wbuf_[0], n_words,
                                &qcache_[0], n_q, 
                                extra + extra_local ) ) {
                        break;
                    }
                    else
                    {
                        for( std::vector< TQNum >::const_iterator i(
                                    marked_list.begin() );
                                i != marked_list.end(); ++i ) {
                            this->queries_.SetMark( *i, true );
                            extra += extra_local;
                        }
                    }
    
                    ++expanded;
                }

                ++curr_idx;
            }

            end_idx_ = curr_idx;
            if( expanded == 0 ) this->end_pass_ = true; // no more subpasses needed
        }

        M_TRACE( CTracer::INFO_LVL, "expanded " << expanded << " queries " );

        // now read the query data once again, loading any extra queries
        // that were marked in the code above
        //
        if( !this->end_pass_ ) {
            curr_idx = 0;
            CReadBinFile in( 
                    this->tmp_store_.Register( CQueryStore::QDUMP_NAME ) );

            while( in.Read( (char *)&q, sizeof( CQueryData ), true ) ) {
                size_t n_words( 3 + ((q.Len() - 1)>>WSHIFT) );
                n_words <<= 1;

                if( in.Read( 
                            (char *)&wbuf_[0], 
                            n_words*sizeof( TWord ), true ) == 0 ) {
                    SRPRISM_ASSERT( false );
                }

                // only add the matked queries that are outside the
                // current subpass index range; the rest have been
                // loaded already; these queries are ignored for the
                // purposes of the main search but are probably needed
                // for possible in-place alignment
                //
                if( curr_idx < start_idx_ || curr_idx >= end_idx_ ) {
                    TQNum qn( q.QNum() );

                    if( this->queries_.IsMarked( qn ) ) {
                        q.SetIgnored( true );

                        if( !this->queries_.AddQueryData( 
                                    &wbuf_[0], n_words, &q, 1 ) ) {
                            SRPRISM_ASSERT( false );
                        }
                    }
                }
    
                ++curr_idx;
            }
        }
    }
    else {
        CQueryData q( 0, 0, 0 );    // query data read from temp file
        size_t curr_idx( 0 ),       // index of the query currently being processed
               expanded( 0 );       // total number of expanded queries in this

        {
            CReadBinFile in( 
                    this->tmp_store_.Register( CQueryStore::QDUMP_NAME ) );
                    
            while( in.Read( (char *)&q, sizeof( CQueryData ), true ) ) {
                size_t n_words( 3 + ((q.Len() - 1)>>WSHIFT) );
                n_words <<= 1;

                // reading raw query data
                //
                if( in.Read( 
                            (char *)&wbuf_[0], 
                            n_words*sizeof( TWord ), true ) == 0 ) {
                    SRPRISM_ASSERT( false );
                }

                // expand query if it belongs to the current subpass
                //
                if( curr_idx >= start_idx_ ) {
                    size_t n_q( 0 ); // number of query data entries resulting
                                     //     from expansion of one query
                    q.SetRawData( &wbuf_[0] );

                    // skip queries that should not be searched in this pass
                    //
                    if( q.IsSALong() ||
                            TBase::template IgnoreForPass< false >( 
                                q, &this->queries_ ) || 
                            this->queries_.template MaxErr< 
                                TScoring, false >( q.QNum() ) <
                                    MinErrForBlowUp( q ) ||
                            this->n_err_ < MinErrForBlowUp( q ) ) {
                        ++curr_idx; continue;
                    }

                    qcache_.clear();
                    BlowUp( q, n_q );
                            
                    // try to append expanded query to the query store
                    //
                    if( !this->queries_.AddQueryData( 
                                &wbuf_[0], n_words, 
                                (qcache_.empty() ? nullptr : &qcache_[0]), n_q ) ) {
                        break;
                    }

                    ++expanded;
                }
    
                ++curr_idx;
            }

            end_idx_ = curr_idx;
            if( expanded == 0 ) this->end_pass_ = true; // no more sub-passes 
                                                        // are needed
        }

        M_TRACE( CTracer::INFO_LVL, "expanded " << expanded << " queries " );
    }
}

//------------------------------------------------------------------------------
//
// create a result structure from the alignment, update housekeeping information
// for the query, and save the result, if necessary
//
template< int search_mode, int hash >
inline void CSearchPass_ByPaired< search_mode, hash, false >::PostProcessMatch( 
        const CHit & hit, int )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;

    if( !this->paired_search_ ) { 
        TPos p( hit.Anchor() );
        TSeqSize sl( hit.AnchorSubjLen() );

        if( hit.Strand() == STRAND_RV ) {
            p -= sl;
        }

        if( !this->seqstore_.CheckRegion( p, sl ) ) return;
    }

    const CQueryData & q( hit.GetQueryData() );
    int min_err( this->queries_.MinErr( q ) );
    TQNum qnum( q.QNum() );

    if( this->queries_.template AddSingleResult< TScoring >( 
                qnum, min_err, hit.AlignLen(), hit.FullNErr(), 
                hit.FullNId(), hit.FullNDel(), hit.FullNGOpen() ) ) {
        TPos encoded_pos( hit.Anchor() );
        CSeqStore::TDecSeqData s_pos( 
                this->seqstore_.DecodePos( encoded_pos ) );
        const CQueryData & q( hit.GetQueryData() );
        TQNum qnum( q.QNum() );
        CResult r( this->tmp_res_mgr_.Save( CResult::EstimateLen(
                        1, hit.FullNErr() ) ) );
        r.Init( 
                qnum, s_pos.first, s_pos.second, hit.AlignLen(),
                hit.LeftOffset(), hit.RightOffset(), 
                hit.Strand(), hit.FullNErr() );
        CHit::TErrors err_info( hit.GetErrors() );

        if( hit.Strand() == seq::STRAND_FW ) {
            for( int i( 0 ); i < hit.FullNErr(); ++i ) {
                r.SetErrorInfo( 
                        0, i, err_info[i].qpos, err_info[i].err_type );
            }
        }
        else {
            int n_err( hit.FullNErr() );

            for( int i( 0 ); i < n_err; ++i ) {
                r.SetErrorInfo( 
                        0, n_err - i - 1, 
                        err_info[i].qpos, err_info[i].err_type );
            }
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash >
inline bool CSearchPass_ByPaired< search_mode, hash, true >::NeedsSearch( 
        const CQueryData & q, TQNum qn ) const
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    return !this->queries_.IsSlave( qn ) && 
                this->queries_.template MaxErr< TScoring, true >( qn ) >= 
                    q.NErr();
}

//------------------------------------------------------------------------------
template< int search_mode, int hash >
inline void CSearchPass_ByPaired< search_mode, hash, true >::SearchInPlace( 
        const CHit & h, TQNum * qs, TQNum * qe, 
        TPos p, int n_err, TAligner & aligner )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    typedef CSeqStore::TDecSeqData TSData;
    int min_err( this->queries_.MinErr( h.GetQueryData() ) );
    bool left_primary_hit( this->queries_.IsLeft( *qs ) );
    size_t ipam_idx( left_primary_hit ? 0 : 2 );
    ipam_idx += (h.Strand() == STRAND_FW) ? 0 : 1;
    T_IPAM ipam( ipam_vec_.data[ipam_idx] );
    TQNum pqn( this->queries_.GetMate( *qs ) );
    if( this->queries_.IsIgnored( pqn ) ) return;
    CQueryData & q( this->queries_.PrimaryData( pqn ) );
    AdjustInPlaceSegs( 
            q.Len(), aligner.GetAnchorStart(), aligner.GetAnchorEnd() );

    // actual requested number of errors depends on whether the query 
    // is < 32 bpp in length and on the current rank of the mate
    //
    {
        int me = MaxQueryErrors( q.IsShort() );
        n_err = std::min( n_err, me );
        int re( this->queries_.template MaxErr< TScoring, true >( pqn ) );
        n_err = std::min( re, n_err );
    }

    // actual in-place search happens here
    //
    results_.clear();
    aligner.Search< SSearchModeTraits< search_mode >::PARTIAL_ALIGNMENT >( 
            q, this->queries_.GetSAStart(), this->queries_.GetSAEnd(),
            n_err, ipam, results_ );

    // save the results if needed
    //
    if( !results_.empty() ) {
        TPos ep( h.Anchor() );
        TSData sep( this->seqstore_.DecodePos( ep ) );
        CHit::TErrors ei( h.GetErrors() );

        for( TAligner::TResults::iterator i = results_.begin();
                i != results_.end(); ++i ) {
            TPos epp( i->Anchor() );
            TSData sepp( this->seqstore_.DecodePos( epp ) );
            CHit::TErrors pei( i->GetErrors() );

            {
                TSeqSize epos( sep.second ),
                         ipos( sepp.second );

                if( h.Strand() == STRAND_RV )
                {
                    epos -= h.AnchorSubjLen();
                    ep -= h.AnchorSubjLen();
                }

                if( i->Strand() == STRAND_RV )
                {
                    ipos -= i->AnchorSubjLen();
                    epp -= i->AnchorSubjLen();
                }

                bool good_template(
                        this->queries_.CheckTemplateConstraints(
                            epos, ipos,
                            h.AnchorSubjLen(), i->AnchorSubjLen() ) );

                if( !good_template ||
                    !this->seqstore_.CheckRegionPair(
                        ep, epp, h.AnchorSubjLen(), i->AnchorSubjLen() ) ) {
                    continue;
                }
            }

            /*
            bool good_template(
                    left_primary_hit ?
                    this->queries_.CheckTemplateConstraints(
                        ep, epp, h.AnchorSubjLen(), i->AnchorSubjLen() ) :
                    this->queries_.CheckTemplateConstraints(
                        epp, ep, i->AnchorSubjLen(), h.AnchorSubjLen() ) );

            if( !good_template ||
                !this->seqstore_.CheckRegionPair( 
                        ep, epp, h.AnchorSubjLen(), i->AnchorSubjLen() ) ) {
                continue;
            }
            */

            /*
            if( !this->seqstore_.CheckRegionPair(
                        ep, epp, h.AnchorSubjLen(), i->AnchorSubjLen() ) ) {
                continue;
            }
            */

            if( this->queries_.template AddPairedResult< TScoring >(
                        qs, qe, min_err,
                        h.AlignLen(), h.FullNErr(), 
                        h.FullNId(), h.FullNDel(), h.FullNGOpen(),
                        i->AlignLen(), i->FullNErr(), 
                        i->FullNId(), i->FullNDel(), i->FullNGOpen() ) ) {
                for( TQNum * qi( qs ); qi != qe; ++qi ) {
                    CResult r( this->tmp_res_mgr_.Save( CResult::EstimateLen(
                                    2, h.FullNErr(), i->FullNErr() ) ) );

                    if( left_primary_hit ) { // reference query is the left one
                        r.Init( 
                                *qi, sep.first, sep.second, sepp.second,
                                h.AlignLen(), i->AlignLen(),
                                h.LeftOffset(), i->LeftOffset(),
                                h.RightOffset(), i->RightOffset(),
                                h.Strand(), i->Strand(),
                                h.FullNErr(), i->FullNErr() );

                        if( h.Strand() == seq::STRAND_FW ) {
                            for( int i( 0 ); i < h.FullNErr(); ++i ) {
                                r.SetErrorInfo( 
                                        0, i, 
                                        ei[i].qpos,
                                        ei[i].err_type );
                            }
                        }
                        else {
                            int n_err( h.FullNErr() );

                            for( int i( 0 ); i < n_err; ++i ) {
                                r.SetErrorInfo( 
                                        0, n_err - i - 1, 
                                        ei[i].qpos,
                                        ei[i].err_type );
                            }
                        }

                        if( i->Strand() == seq::STRAND_FW ) {
                            for( int j( 0 ); j < i->FullNErr(); ++j ) {
                                r.SetErrorInfo( 
                                        1, j, 
                                        pei[j].qpos,
                                        pei[j].err_type );
                            }
                        }
                        else {
                            int n_err( i->FullNErr() );

                            for( int i( 0 ); i < n_err; ++i ) {
                                r.SetErrorInfo( 
                                        1, n_err - i - 1, 
                                        pei[i].qpos,
                                        pei[i].err_type );
                            }
                        }
                    }
                    else {
                        r.Init(
                                this->queries_.GetMate( *qi ), sep.first,
                                sepp.second, sep.second,
                                i->AlignLen(), h.AlignLen(),
                                i->LeftOffset(), h.LeftOffset(),
                                i->RightOffset(), h.RightOffset(),
                                i->Strand(), h.Strand(),
                                i->FullNErr(), h.FullNErr() );

                        if( i->Strand() == seq::STRAND_FW ) {
                            for( int j( 0 ); j < i->FullNErr(); ++j ) {
                                r.SetErrorInfo( 
                                        0, j, 
                                        pei[j].qpos,
                                        pei[j].err_type );
                            }
                        }
                        else {
                            int n_err( i->FullNErr() );

                            for( int i( 0 ); i < n_err; ++i ) {
                                r.SetErrorInfo( 
                                        0, n_err - i - 1, 
                                        pei[i].qpos,
                                        pei[i].err_type );
                            }
                        }

                        if( h.Strand() == seq::STRAND_FW ) {
                            for( int i( 0 ); i < h.FullNErr(); ++i ) {
                                r.SetErrorInfo( 
                                        1, i, 
                                        ei[i].qpos,
                                        ei[i].err_type );
                            }
                        }
                        else {
                            int n_err( h.FullNErr() );

                            for( int i( 0 ); i < n_err; ++i ) {
                                r.SetErrorInfo( 
                                        1, n_err - i - 1, 
                                        ei[i].qpos,
                                        ei[i].err_type );
                            }
                        }
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
//
// performs the in-place alignment for the mates of the query and
// its duplicates
//
template< int search_mode, int hash >
inline void CSearchPass_ByPaired< search_mode, hash, true >::PostProcessMatch( 
        const CHit & hit, int n_err )
{
    typedef CSeqStore::TDecSeqData TSData;
    const CQueryData & q( hit.GetQueryData() );
    TQNum qn( q.QNum() );
    bool skip( false );
    TSegDescr segs[2];
    TPos p( hit.Anchor() );

    if( this->queries_.IsUnique( qn ) && !NeedsSearch( q, qn ) ) {
        skip = true;
    }

    // generate 2 search regions
    {
        TSeqSize seglen( 2*this->queries_.PairFuzz() ),
                 segoff( this->queries_.PairDistance() - 
                            this->queries_.PairFuzz() );
        if( hit.Strand() == seq::STRAND_RV ) p -= hit.AnchorSubjLen();
        TSData subj_data( this->seqstore_.DecodePos( p ) );

        // left region
        segs[0].first = segs[0].second = 0;
        TSeqSize ee( subj_data.second + hit.AnchorSubjLen() ),
                 end( ee - std::min( segoff, ee ) );

        if( end > 0 ) {
            TSeqSize l( std::min( end, seglen ) );
            segs[0].first = this->seqstore_.EncodePos( 
                    std::make_pair( subj_data.first, end - l ) );
            segs[0].second = std::min( ee - segs[0].first, l + n_err );
        }

        // right region
        segs[1].second = 0;
        TSeqSize tail( this->seqstore_.FwTailLen( p ) );

        if( tail >= segoff ) {
            segs[1].first = p + segoff;
            segs[1].second = std::min( seglen, tail - segoff );
            segs[1].first -= std::min( (TSeqSize)n_err, segoff );
            segs[1].second += std::min( (TSeqSize)n_err, segoff );
        }
    }

    if( !skip ) {
        //
        // attempt in-place search
        //
        TAligner aligner( 
                this->seqstore_, ip_segs_, p, p + hit.AnchorSubjLen(), this->ip_ma_, 
                this->queries_.ResLimit(),
                &this->pass_stats_.n_aligns, 
                &this->pass_stats_.n_ualigns,
                &this->pass_stats_.n_inplace,
                &this->pass_stats_.n_inplace_aligns );

        if( this->queries_.IsUnique( qn ) ) {
            // generate initial search intervals
            //
            ip_segs_.clear();

            if( segs[0].second > 0 ) {
                ip_segs_.push_back( 
                        TAligner::Seg( segs[0].first, segs[0].second, true ) );
            }

            if( segs[1].second > 0 ) {
                ip_segs_.push_back( 
                        TAligner::Seg( segs[1].first, segs[1].second, false ) );
            }

            SearchInPlace( hit, &qn, &qn + 1, p, n_err, aligner );
        }
        else {
            TQNum * dqs( this->queries_.DupStart( qn ) ),
                  * dqe( this->queries_.DupEnd( qn ) );

            for( TQNum * dqi( dqs ); dqi != dqe; ) {
                if( !this->queries_.Done4Search( *dqi ) && 
                        NeedsSearch( q, *dqi ) ) {
                    TQNum * dqie( dqi );

                    // the set of pairs where both mates are the same
                    // is searched only once
                    //
                    while( dqie != dqe && 
                            this->queries_.EquivMates( *dqi, *dqie ) ) {
                        ++dqie;
                    }

                    // generate initial search intervals
                    //
                    ip_segs_.clear();

                    if( segs[0].second > 0 ) {
                        ip_segs_.push_back( TAligner::Seg( 
                                    segs[0].first, segs[0].second, true ) );
                    }

                    if( segs[1].second > 0 ) {
                        ip_segs_.push_back( TAligner::Seg( 
                                    segs[1].first, segs[1].second, false ) );
                    }

                    SearchInPlace( hit, dqi, dqie, p, n_err, aligner );
                    dqi = dqie;
                }
                else ++dqi;
            }
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired >
inline void CSearchPass< search_mode, hash, paired >::ProcessCandidate( 
        const CQueryData & q, TPos pos, int n_err, TStrand s )
{
    ++this->pass_stats_.n_candidates;

    // apply a hash-specific aligner to the (q,pos) and if successful
    // apply a post-aligner; post-aligner for unpaired reads just saves
    // the alignment, while post-aligner for paired reads first attempts
    // an in-place alignment of the mate and depending on the outcome
    // saves either paired or unpaired result
    //
    TSeqSize qoff( q.GetOff() );
    CHit hit( q, s );
    this->template GenBreakSegs< hash >( q );

    if( hit.Extend< SSearchModeTraits< search_mode >::PARTIAL_ALIGNMENT >( 
                this->seqstore_, this->main_ma_, 
                this->queries_.GetSAStart(), this->queries_.GetSAEnd(),
                qoff, pos, HASH_LEN, n_err, &this->bsegs_, 
                &this->pass_stats_.n_aligns,
                &this->pass_stats_.n_ualigns ) ) {
        TBaseByPaired::PostProcessMatch( hit, n_err );
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired >
inline void CSearchPass< search_mode, hash, paired >::ProcessQuery( 
        const CQueryData & query )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    typedef CIndexIterator::TPosIter TIter;
    TStrand qstrand = query.Strand();

    // adjust the max number of errors for the search: short queries
    // have different target number of errors; also if we already found
    // better results than 'n_err_' reduce target number of errors
    // accordingly
    //
    int n_err = std::min( this->n_err_, MaxQueryErrors( query.IsShort() ) );
    n_err = std::min( 
            n_err, 
            this->queries_.template GroupMaxErr< TScoring, paired >( 
                query.QNum() ) );

    // process subject position lists for each strand
    //
    if( this->randomize_ ) {
        TIter s1( this->idx_->PosStart( STRAND_FW ) ),
              e1( this->idx_->PosEnd( STRAND_FW ) ),
              s2( this->idx_->PosStart( STRAND_RV ) ),
              e2( this->idx_->PosEnd( STRAND_RV ) );
        Uint4 sz1( e1 - s1 ), sz( e2 - s2 + sz1 );

        if( sz > 0 ) {
            std::vector< Uint4 > rnd_map_;
            rnd_map_.reserve( sz );
            for( size_t i( 0 ); i < sz; ++i ) rnd_map_.push_back( i );

            for( size_t i( 0 ); i < sz - 1; ++i ) {
                // Uint4 idx( random()%(sz - i) );
                Uint4 idx( rand()%(sz - i) );
                std::swap( rnd_map_[idx], rnd_map_[sz - i - 1] );
            }

            for( size_t i( 0 ); i < sz; ++i ) {
                Uint4 idx( rnd_map_[i] );

                if( idx < sz1 ) {
                    TStrand s( CombineStrands( qstrand, STRAND_FW ) );
                    ProcessCandidate( query, *(s1 + idx), n_err, s );
                }
                else {
                    TStrand s( CombineStrands( qstrand, STRAND_RV ) );
                    ProcessCandidate( query, *(s2 + (idx - sz1)), n_err, s );
                }
            }
        }
    }
    else {
        for( TStrand j = 0; j < N_STRANDS; ++j ) {
            TStrand strand( CombineStrands( qstrand, j ) );

            for( TIter i = this->idx_->PosStart( j ); 
                    i != this->idx_->PosEnd( j ); ++i ) {
                ProcessCandidate( query, *i, n_err, strand );
            }
        }
    }

    if( hash == HASH_NORMAL ) {
        this->queries_.SetHashDone( query.QNum(), query.GetCurrHashIdx() );
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired >
void CSearchPass< search_mode, hash, paired >::ProcessQueryBlock( 
        const CQueryData & qstart, const CQueryData & qend )
{
    const CQueryData * start( &qstart ), * end( &qend );

    for( const CQueryData * i( start ); i != end; ++i ) {
        if( !this->queries_.template GroupDone4Search< paired >( 
                    i->QNum() ) ) {
            ProcessQuery( *i );
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired > inline void 
CSearchPass< search_mode, hash, paired >::ProcessSpecialQueryBlockPos( 
        CQueryData * qbase, 
        typename TBase::TQExtDataTable & edt, 
        const typename TBase::TBNF::TQExtSet & q_ext_set, 
        TPos pos, int n_err, TStrand s )
{
    typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
    typedef typename TBase::TBNF::TQExtSet::const_iterator TQExtIter;

    // search queries
    //
    for( TQExtIter qei( q_ext_set.begin() ); qei != q_ext_set.end(); ++qei ) {
        const typename TBase::SQExtData & qed( edt[*qei] );
        const CQueryData * qs( qbase + qed.idx ), * qe( qs + qed.count );

        for( ; qs != qe; ++qs ) {
            if( !this->queries_.template GroupDone4Search< paired >( 
                        qs->QNum() ) ) {
                int n_err_lcl( std::min( 
                            n_err,
                            this->queries_.template GroupMaxErr< 
                                TScoring, paired >( qs->QNum() ) ) );
                ProcessCandidate( *qs, pos, n_err_lcl, s );
            }
        }
    }

    // purge queries
    //
    int filter_n_err( qbase->NErr() + this->bnf_.GetNErr() );

    for( TQExtIter qei( q_ext_set.begin() ); qei != q_ext_set.end(); ++qei ) {
        typename TBase::SQExtData & qed( edt[*qei] );
        CQueryData * qs( qbase + qed.idx ), * qi( qs );

        while( qi - qs < qed.count ) {
            TQNum qn( qi->QNum() );

            if( this->queries_.template GroupMaxErr< 
                        TScoring, paired >( qn ) < filter_n_err ||
                    this->queries_.template GroupDone4Search< paired >( qn ) ) {
                --qed.count;
                this->queries_.SwapQueryData( qi, qs + qed.count );
            }
            else ++qi;
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired > inline void 
CSearchPass< search_mode, hash, paired >::ProcessSpecialQueryBlock_ExtExt( 
        CQueryData * qbase, const CIndexIterator::TExtData & ext, 
        typename TBase::TQExtDataTable & edt, 
        int n_err, int n_err_filter, TSeqSize ext_len, 
        TStrand qstrand, TStrand estrand )
{
    typedef typename TBase::TBNF::TMatchSet TMatchSet;
    typedef CIndexIterator::TPosIter TPosIter;

    this->bnf_.Init( 
            &ext[0], &edt[0], ext.size(), edt.size(), ext_len, n_err_filter );

    std::vector< Uint4 > rnd_map_;

    for( int cnerr( 0 ); cnerr <= n_err_filter; ++cnerr ) {
        while( this->bnf_.Next() >= 0 ) {
            const TMatchSet & ms( this->bnf_.GetMatchSet() );
            TPosIter pos( this->idx_->PosStart( estrand ) + 
                            ext[ms.s_ext].index ),
                     epos( pos + ext[ms.s_ext].fnpos );

            if( this->randomize_ ) {
                TPosIter rpos( epos );
                epos += ext[ms.s_ext].rnpos;
                if( epos == pos ) continue;
                size_t sz( epos - pos ), sz1( rpos - pos );
                rnd_map_.clear();
                rnd_map_.reserve( sz);
                for( size_t i( 0 ); i < sz; ++i ) rnd_map_.push_back( i );

                for( size_t i( 0 ); i < sz - 1; ++i ) {
                    // Uint4 idx( random()%(sz - i) );
                    Uint4 idx( rand()%(sz - i) );
                    std::swap( rnd_map_[idx], rnd_map_[sz - i - 1] );
                }

                for( size_t i( 0 ); i < sz; ++i ) {
                    Uint4 idx( rnd_map_[i] );
                    TStrand s( CombineStrands( 
                                qstrand, (idx < sz1) ? STRAND_FW 
                                                     : STRAND_RV ) );
                    TPosIter p( pos + rnd_map_[i] );
                    ProcessSpecialQueryBlockPos( 
                            qbase, edt, ms.q_ext_set, *p, n_err, s );
                }
            }
            else {
                TStrand s( CombineStrands( qstrand, STRAND_FW ) );

                for( ; pos != epos; ++pos ) {
                    ProcessSpecialQueryBlockPos( 
                            qbase, edt, ms.q_ext_set, *pos, n_err, s );
                }

                epos += ext[ms.s_ext].rnpos;
                s = CombineStrands( qstrand, STRAND_RV );

                for( ; pos != epos; ++pos ) {
                    ProcessSpecialQueryBlockPos( 
                            qbase, edt, ms.q_ext_set, *pos, n_err, s );
                }
            }
        }

        // purge queries
        //
        typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;
        int pnerr( qbase->NErr() + this->bnf_.GetNErr() );

        for( size_t j( 0 ); j < edt.size(); ++j ) {
            typename TBase::SQExtData & qed( edt[j] );
            CQueryData * qs( qbase + qed.idx ), * qi( qs );

            while( qi - qs < qed.count ) {
                TQNum qn( qi->QNum() );

                if( this->queries_.template GroupMaxErr< 
                            TScoring, paired >( qn ) < pnerr ||
                        this->queries_.template GroupDone4Search< 
                            paired >( qn ) ) {
                    --qed.count;
                    this->queries_.SwapQueryData( qi, qs + qed.count );
                }
                else ++qi;
            }
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired > inline void 
CSearchPass< search_mode, hash, paired >::ProcessSpecialQueryBlock_Ext( 
        CQueryData & qstart, CQueryData & qend, 
        int n_err, int n_err_filter, TSeqSize ext_len,
        TStrand qstrand, TStrand estrand )
{
    const CIndexIterator::TExtData & extensions(
            this->idx_->Extensions( estrand ) );
    if( extensions.size() == 0 ) return;
    CQueryData * qbase( &qstart ), * qendptr( &qend );

    typedef typename TBase::SQExtData TQExtData;
    typename TBase::TQExtDataTable & edt( this->ext_data_table_ );

    while( qbase != qendptr ) {
        //
        // n_q is the size of block of hashes to process
        //
        size_t n_q( std::min( 
                    (size_t)(qendptr - qbase), 
                    TBase::MAX_QEXT_DATA_IDX + 1 ) );
        CQueryData * qi( qbase ), * qe( qbase + n_q );

        while( qi != qe ) {
            edt.clear();
            
            while( qi != qe && (size_t)(qi - qbase) < TBase::QEXT_TABLE_SIZE ) {
                //
                // initialize and process a block of uniq'ed extensions
                //
                TWord ext( qi->GetExtension() );
                edt.push_back( TQExtData( ext, qi - qbase ) );
                TQExtData & d( *edt.rbegin() );
                ++qi;

                while( qi != qe &&
                        (size_t)(qi - qbase) < TBase::QEXT_TABLE_SIZE &&
                        ext == qi->GetExtension() ) {
                    ++d.count;
                    ++qi;
                }
            }
                
            ProcessSpecialQueryBlock_ExtExt( 
                    qbase, extensions, edt, n_err, n_err_filter, 
                    ext_len, qstrand, estrand );
            qbase = qi;
        }
    }
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired >
void CSearchPass< search_mode, hash, paired >::ProcessSpecialQueryBlock( 
        CQueryData & qstart, CQueryData & qend )
{
    //
    // this function tries to reduce the complexity of search for very
    // repetitive n-mers by first trying to match their 16-base extensions
    // and only considering pos/query pairs for which extesnions match
    // for further search
    //

    int n_err( std::min( this->n_err_, MaxQueryErrors( qstart.IsShort() ) ) );

    // TODO: Check if this is really needed, i.e. it is possible that
    //       query blowup for short queries never results in 2 errors
    //       in the seed.
    //
    if( n_err >= qstart.NErr() ) {
        // number of errors for matching subject and query extensions
        //
        int n_err_filter = std::min( 2, n_err - qstart.NErr() );

        // extension length to work with: normally 16, but can be less for
        // short queries
        //
        TSeqSize ext_len( qstart.ExtLen() );

        if( qstart.IsPalindrome() ) {
            if( qstart.IsRightExtDir() ) {
                ProcessSpecialQueryBlock_Ext( 
                        qstart, qend, n_err, n_err_filter, ext_len, 
                        STRAND_FW, STRAND_FW );
                ProcessSpecialQueryBlock_Ext( 
                        qstart, qend, n_err, n_err_filter, ext_len, 
                        STRAND_RV, STRAND_RV );
            }
            else {
                ProcessSpecialQueryBlock_Ext( 
                        qstart, qend, n_err, n_err_filter, ext_len, 
                        STRAND_FW, STRAND_RV );
                ProcessSpecialQueryBlock_Ext( 
                        qstart, qend, n_err, n_err_filter, ext_len, 
                        STRAND_RV, STRAND_FW );
            }
        }
        else {
            TStrand qstrand( qstart.Strand() ), 
                    estrand( qstart.IsRightExtDir() ? 
                                qstrand : ReverseStrand( qstrand ) );
            ProcessSpecialQueryBlock_Ext(
                    qstart, qend, n_err, n_err_filter, ext_len,
                    qstrand, estrand );
        }
    }

    TBase::template SetHashDone< hash >( qstart, qend );
}

//------------------------------------------------------------------------------
template< int search_mode, int hash, bool paired >
void CSearchPass< search_mode, hash, paired >::Run(void)
{
    while( !this->end_pass_ ) {
        this->pass_stats_.Clean();
        this->queries_.StartProcess();
        CQueryStore::CEquivIterator eq_iter( this->queries_.EquivIterator() );

        {
            M_TRACE( CTracer::INFO_LVL, eq_iter.Total() << " primary queries" );
            this->idx_.reset( new CIndexIterator( this->idx_basename_ ) );
            std::ostringstream os;
            bool idx_done( false );

            while( eq_iter.Next() ) {
                CQueryData & qstart( eq_iter.Start() );
                TPrefix qprefix = qstart.Prefix();
                CQueryData & qend( eq_iter.End() );

                if( !idx_done && !this->idx_->Seek( qprefix ) ) { 
                    idx_done = true;
                }

                if( idx_done ) {
                    TBase::template SetHashDone< hash >( qstart, qend );
                    continue;
                }

                if( (TWord)this->idx_->Prefix() == qprefix ) {
                    if( this->repeat_threshold_ == 0 ||
                            this->idx_->NPos() <= this->repeat_threshold_ ) {
                        if( this->idx_->Special() ) { 
                            if( !skip_bad_ ) {
                                ProcessSpecialQueryBlock( qstart, qend ); 
                            }

                            TBase::SetBadHash( qstart, qend );
                        }
                        else if( !skip_good_ ) {
                            ProcessQueryBlock( qstart, qend );
                        }
                    }
                    else {
                        if( !paired ) {
                            for( CQueryData * qi( &qstart ); 
                                    qi != &qend; ++qi ) {
                                if( hash == HASH_NORMAL ) {
                                    this->queries_.IncrRepHashCount( 
                                            qi->QNum() );
                                }
                                else this->queries_.SetRepBUHash( qi->QNum() );
                            }
                        }
                        else if( hash != HASH_NORMAL ) {
                            for( CQueryData * qi( &qstart ); 
                                    qi != &qend; ++qi ) {
                                this->queries_.SetRepBUHash( qi->QNum() );
                            }
                        }
                    }
                }
                else TBase::template SetHashDone< hash >( qstart, qend );
            }
        }

        this->pass_stats_.UpdateGlobalStats();
        this->pass_stats_.Report();
        this->NextSubPass();
    }
}

END_NS( srprism )
END_STD_SCOPES

