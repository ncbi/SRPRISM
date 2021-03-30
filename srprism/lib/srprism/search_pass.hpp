/*  $Id: search_pass.hpp 477203 2015-08-27 13:22:18Z morgulis $
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

#ifndef __SRPRISM_SEARCH_PASS_HPP__
#define __SRPRISM_SEARCH_PASS_HPP__

#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <memory>

#include "../common/exception.hpp"
#include "../common/tmpstore.hpp"
#include "stat.hpp"
#include "memmgr.hpp"
#include "seqstore.hpp"
#include "tmpres_mgr.hpp"
#include "index_iterator.hpp"
#include "align.hpp"
#include "scratch.hpp"
#include "query_data.hpp"
#include "query_store.hpp"
#include "seqiter.hpp"
#include "inplace_align.hpp"
#include "search_mode.hpp"
#include "bnf.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <memory>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/tmpstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/tmpres_mgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/index_iterator.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/align.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/scratch.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_data.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_store.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqiter.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/inplace_align.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/search_mode.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/bnf.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

class CSearchPassDef
{
    public:

        // initalizer of a search pass
        struct SInitData
        {
            std::string index_basename;     // base name of the index data
            common::Uint4 res_limit;        // limit on the number of reported 
                                            //      results
            common::Uint4 repeat_threshold; // do not process 16-mers with more
                                            //      than this many occurences
                                            //      in the database
            common::Uint2 pair_distance;    // target distance between mates
                                            //      for paired searches
            common::Uint2 pair_fuzz;        // how much deviation from the
                                            //      target distance is allowed
            common::Uint1 n_err;            // search for alignments with at
                                            //      most this many errors

            S_IPAM ipam_vec;                // paired result configuration
                                            // spec

            CMemoryManager * mem_mgr_p;         // object to track memory usage
            CTmpResMgr * tmpres_mgr_p;          // temporary storage for results
            CSeqStore * seqstore_p;             // subject sequence data storage
            common::CTmpStore * tmp_store_p;    // temporary file name manager

            CScratchBitMap * scratch_p; // scratch space to record visited
                                        //      subject-query extension pairs
                                        //      in bad 16-mer cases
            CQueryStore * queries_p;    // query data manager
            CStatMap * search_stats;    // global search statistics
            bool paired_search;         // indication of whether search as a whole
                                        // is on paired queries
            bool randomize;             // randomize results on subject pos
            bool random_seed;           // use random seed to initialize RNG.
        };

};

//------------------------------------------------------------------------------
// most common CSearchPass functionality
//
template< int search_mode >
class CSearchPass_Base
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode PASS_SKIP = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == PASS_SKIP ) return "skipping pass";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

    protected:

        typedef typename SSearchModeTraits< search_mode >::TScoringSys TScoring;

        // aliases
        typedef CIndexIterator::TPos TPos;

        // name of the temporary file containing the queries of the current
        // batch in binary form
        static const char * PQDUMP_NAME;

        // Instance constructor
        //
        CSearchPass_Base( const CSearchPassDef::SInitData & init_data );

        // helper method to update hash use information for a range
        // of queries
        //
        template< int hash >
        void SetHashDone( const CQueryData & qstart, const CQueryData & qend )
        {
            if( hash == HASH_NORMAL ) {
                const CQueryData * end( &qend );

                for( const CQueryData * i( &qstart ); i != end; ++i ) {
                    queries_.SetHashDone( i->QNum(), i->GetCurrHashIdx() );
                }
            }
        }

        static void SetBadHash( CQueryData & qstart, CQueryData & end )
        { for( CQueryData * i( &qstart ); i != &end; ++i ) i->SetBadHash(); }

        //
        // fill extension break segments by hash and query
        //
        template< int hash_type > void GenBreakSegs( const CQueryData & );

        //######################################################################

        //
        // check whether a query should be ignored for a particular hash
        //
        // ignored if it is marked as 'done' in the query store or if
        // results with high rank are already found for it
        //
        template< bool paired > static bool 
        IgnoreForPass( const CQueryData & qh, const CQueryStore * queries_p )
        {
            return queries_p->GroupDone4Search< paired >( qh.QNum() ) ||
                (queries_p->template GroupMaxErr< TScoring, paired >( 
                    qh.QNum() ) < queries_p->MinErr( qh ));
        }

        CMemoryManager & mem_mgr_;  // memory manager
        CTmpResMgr & tmp_res_mgr_;  // temporary result storage
        CSeqStore & seqstore_;      // subject sequence data storage
        CQueryStore & queries_;     // query data manager
        CTmpStore & tmp_store_;     // temporary file name manager

        CExtensionSpaceAllocator main_ma_;  // matrix allocator for initial alignments
        CExtensionSpaceAllocator ip_ma_;    // matrix allocator for inplace alignments

        std::auto_ptr< CIndexIterator > idx_; // index iterator
        std::string idx_basename_;            // index base name

        size_t res_limit_;          // limit on the number of reported results 
                                    //      per query
        size_t repeat_threshold_;   // do not consider hash values with more
                                    //      than this many occurences in the 
                                    //      database
        TSeqSize pair_distance_;    // target distance between mates
        TSeqSize pair_fuzz_;        // allowed deviation from 'pair_distance_'
        int n_err_;                 // target max number of errors
        bool end_pass_;             // flag indicating end of pass 
                                    //      (no more sub-passes)
        bool paired_search_;        // true if the search as a whole is 
                                    //      performed on paired queries
        CBreakSegs bsegs_;          // break segment container for extensions

        // min number of hits guaranteed by a matrix allocator before
        // cleanup
        //
        static const size_t MATRIX_ALLOC_INIT_HITS = 32*1024ULL;

        /////////////////////////////////////////////////////////
        //
        // DECLARATIONS TO SUPPORT BAD 16-MER FILTER CODE
        //
        /////////////////////////////////////////////////////////
        struct SQExtData
        {
            SQExtData( TWord e = 0, common::Uint2 i = 0 )
                : ext( e ), count( 1 ), idx( i )
            {}

            TWord ext;
            common::Uint2 count;
            common::Uint2 idx;
        };

        typedef CBadNMerFilter< CIndexIterator::TExtInfo, SQExtData > TBNF;
        typedef std::vector< SQExtData > TQExtDataTable;

        static const size_t QEXT_TABLE_SIZE = TBNF::MAX_QDATA_SIZE;
        static const size_t MAX_QEXT_DATA_IDX = 
            common::SIntTraits< common::Uint2 >::MAX;

        bool AlignFilter( 
                TExtension e1, TExtension e2, TWord m, 
                TSeqSize ext_len, int n_err )
        {
            if( e1 == e2 && m == 0 ) return false;
            ++pass_stats_.n_filter;

            static const size_t LBITS = 
                seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
            static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

            TExtension mask( 
                    common::MaskFront< TExtension >( ext_len<<LSHIFT ) );
            e1 &= mask; e2 &= mask; mask <<= LBITS;

            if( m == 0 ) {
                if( n_err == 1 ) {
                    return !CFastAlignCheck< TExtension, 1 >()( e1, e2, mask );
                }
                else return !CFastAlignCheck< TExtension, 2 >()( e1, e2, mask );
            }
            else {
                if( n_err == 1 ) {
                    return !CFastAlignCheckWithMasks< TExtension, 1 >()(
                            e1, e2, m, 0, mask );
                }
                else {
                    return !CFastAlignCheckWithMasks< TExtension, 2 >()(
                            e1, e2, m, 0, mask );
                }
            }
        }

        TQExtDataTable ext_data_table_; // query extension working set

        struct _SPassStats
        {
            _SPassStats( CStatMap * search_stats )
            {
                Clean();
                global_n_aligns = search_stats->GetCounter( STAT_N_ALIGNS );
                global_n_ualigns = search_stats->GetCounter( STAT_N_UALIGNS );
                global_n_filter = search_stats->GetCounter( STAT_N_FILTER );
                global_n_candidates = 
                    search_stats->GetCounter( STAT_N_CANDIDATES );
                global_n_inplace = search_stats->GetCounter( STAT_N_INPLACE );
                global_n_inplace_aligns = 
                    search_stats->GetCounter( STAT_N_INPLACE_ALIGNS );
            }

            void Clean( void )
            {
                n_aligns = n_ualigns = n_filter = n_candidates 
                         = n_inplace = n_inplace_aligns = 0;
            }

            void UpdateGlobalStats( void ) 
            {
                *global_n_aligns         += n_aligns;
                *global_n_ualigns        += n_ualigns;
                *global_n_filter         += n_filter;
                *global_n_candidates     += n_candidates;
                *global_n_inplace        += n_inplace;
                *global_n_inplace_aligns += n_inplace_aligns;
            }

            void Report( void )
            {
                M_TRACE( common::CTracer::INFO_LVL, "pass statistics: " );
                M_TRACE( common::CTracer::INFO_LVL, 
                         "\textensions:                " << n_aligns << 
                         " (" << *global_n_aligns  << ")" );
                M_TRACE( common::CTracer::INFO_LVL,
                         "\tunidirectional extensions: " << n_ualigns << 
                         " (" << *global_n_ualigns << ")" );
                M_TRACE( common::CTracer::INFO_LVL, 
                         "\tbad nmer extension checks: " << n_filter << 
                         " (" << *global_n_filter  << ")" );
                M_TRACE( common::CTracer::INFO_LVL, 
                         "\tresult candidates:         " << n_candidates << 
                         " (" << *global_n_candidates  << ")" );
                M_TRACE( common::CTracer::INFO_LVL, 
                         "\tin-place scans:            " << n_inplace << 
                         " (" << *global_n_inplace  << ")" );
                M_TRACE( common::CTracer::INFO_LVL, 
                         "\tin-place seeds:            " << n_inplace_aligns << 
                         " (" << *global_n_inplace_aligns  << ")" );
            }

            size_t * global_n_aligns;
            size_t * global_n_ualigns;
            size_t * global_n_filter;
            size_t * global_n_candidates;
            size_t * global_n_inplace;
            size_t * global_n_inplace_aligns;

            size_t n_aligns;
            size_t n_ualigns;
            size_t n_filter;
            size_t n_candidates;
            size_t n_inplace;
            size_t n_inplace_aligns;
        } pass_stats_;

        TBNF bnf_;
        bool randomize_;
        bool random_seed_;
};

//------------------------------------------------------------------------------
template< int search_mode > template< int hash > inline void
CSearchPass_Base< search_mode >::GenBreakSegs( const CQueryData & q )
{
    if( hash == HASH_NORMAL ) {
        bsegs_.Clear();
        int n_hashes( q.GetNHashes() );
        TQNum qnum( q.QNum() );
        bool bad_check( n_hashes > 2 && q.IsSALong() );

        for( int i( 0 ); i < n_hashes; ++i ) {
            if( queries_.HashDone( qnum, i ) ) {
                TSeqSize off( q.GetHashOffset( i ) );

                if( q.GetBadHash( i ) && bad_check ) {
                    if( q.RightExtDir( i ) ) {
                        bsegs_.Add( off, off + HASH_LEN, 1, bsegs_.Size() + 1 );
                        off += HASH_LEN;
                        bsegs_.Add( off, off + HASH_LEN, 3, bsegs_.Size() - 1 );
                    }
                    else {
                        bsegs_.Add( off - HASH_LEN, off, 3, bsegs_.Size() + 1 );
                        bsegs_.Add( off, off + HASH_LEN, 1, bsegs_.Size() - 1 );
                    }
                }
                else bsegs_.Add( off, off + HASH_LEN, 1 );
            }
        }
    }
    else if( hash == HASH_BLOWUP ) {
        SRPRISM_ASSERT( q.GetNHashes() < 4 && !q.IsSALong() );
        bsegs_.Clear();

        if( q.AmbigSA() ) {
            int bu_hash_idx( q.GetBUHashIdx() );
    
            if( bu_hash_idx != 3 ) {
                int n_hashes( q.GetNHashes() );
            
                for( int i( 0 ); i < n_hashes; ++i ) {
                    if( !q.HashAmbig( i ) ) {
                        TSeqSize s( q.GetHashOffset( i ) );
                        bsegs_.Add( s, s + HASH_LEN, 1 );
                    }
                }
            }
        }
        else if( !q.IsShort() ) {
            TSeqSize s0( q.GetHashOffset( 0 ) ),
                     s1( q.GetHashOffset( 1 ) ),
                     s2( q.GetHashOffset( 2 ) );

            if( (s0 == s1 || s1 == s2) && q.NErr() > 1 ) return;

            if( s1 - s0 == HASH_LEN ) bsegs_.Add( s0, s0 + HASH_LEN, 1 );
            else {
                TSeqSize s2( q.GetHashOffset( 2 ) );
                bsegs_.Add( s2, s2 + HASH_LEN, 1 );
            }
        }
    }
    else SRPRISM_ASSERT( false );
}

//------------------------------------------------------------------------------
// CSearchPass functionality that must by specialized by hash
//
template< int search_mode, int hash, bool paired > class CSearchPass_ByHash;

//##############################################################################
template< int search_mode, bool paired >
class CSearchPass_ByHash< search_mode, HASH_NORMAL, paired > 
    : public virtual CSearchPass_Base< search_mode >
{
    private:

        typedef CSearchPass_Base< search_mode > TBase;

    protected:

        // function object to set up query data structure for a normal hash
        //
        struct HashSetup
        {
            void operator()( 
                    CQueryData & qh, const CQueryStore * queries_p ) const
            {
                int hash_idx( qh.GetCurrHashIdx() );
                TSeqSize hash_start( qh.GetHashOffset( hash_idx ) );
                CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si( 
                        qh.Data(), hash_start, hash_start + HASH_LEN );
                TPrefix prefix( (*si).first ), rc_prefix( 0 );
                ReverseComplement< SEQDATA_CODING >( rc_prefix, prefix );
                qh.SetPrefix( (prefix <= rc_prefix) ? prefix : rc_prefix );
                qh.SetStrand( (prefix < rc_prefix) ? STRAND_FW : STRAND_RV );
                qh.SetPalindrome( prefix == rc_prefix );
                qh.SetRightExtDir( qh.RightExtDir() );
                qh.SetExtLen();

                if( TBase::template IgnoreForPass< paired >( qh, queries_p ) ) {
                    qh.SetIgnored( true );
                }
            }
        };

        // instance constructor
        // 
        CSearchPass_ByHash( const CSearchPassDef::SInitData & init_data )
            : TBase( init_data )
        {
            this->queries_.ClearMarks();
            this->queries_.StartUpdate();
            this->queries_.ForEach( HashSetup() );
            this->queries_.SortQueries();
            this->queries_.SetupCrossLinks();
        }

        // prepare for the next sub-pass 
        // (no additional sub-passes normal hashes)
        //
        void NextSubPass( void ) { this->end_pass_ = true; }
};

//##############################################################################
template< int search_mode, bool paired >
class CSearchPass_ByHash< search_mode, HASH_BLOWUP, paired >
    : public virtual CSearchPass_Base< search_mode >
{
    private:

        typedef CSearchPass_Base< search_mode > TBase;

        // for this hash alignments with at least this many errors 
        // are searched
        //
        static const int MIN_ERR = 1;

    protected:

        // expand the query's 16-base prefix with 1 error
        //
        void BlowUp( const CQueryData & q, size_t & n_q )
        {
            if( q.AmbigSA() ) {
                int hash_idx( q.GetBUHashIdx() );

                if( hash_idx != 3 ) {
                    TSeqSize hash_off( q.GetHashOffset( hash_idx ) ),
                             pos( q.GetAmbigBUPos() );

                    CQueryData::TBlowUpIterator bui( q.GetBlowUpIterator( 
                                hash_off, pos, pos + 1, 
                                q.RightBUDir(), true ) );

                    if( !bui.End() ) {
                        do {
                            CQueryData nq( q );
                            nq.SetPrefix( bui.GetHashValue() );
                            nq.SetType1( bui.GetType_1() );
                            nq.SetType2( bui.GetType_2() );
                            nq.SetPos1( bui.GetPos_1() );
                            nq.SetPos2( bui.GetPos_2() );
                            AddExtra( nq, n_q );
                        }
                        while( bui.Next() );
                    }
                }

                return;
            }

            TSeqSize start( 0 ), end( 0 );
            
            if( q.GetBlowUpParams( start, end ) ) {
                CQueryData::TBlowUpIterator bui( q.GetBlowUpIterator(
                            q.GetHashOffset( 3 ), 
                            start, end, q.RightBUDir(), false ) );

                if( !bui.End() ) {
                    do {
                        CQueryData nq( q );
                        nq.SetPrefix( bui.GetHashValue() );
                        nq.SetType1( bui.GetType_1() );
                        nq.SetType2( bui.GetType_2() );
                        nq.SetPos1( bui.GetPos_1() );
                        nq.SetPos2( bui.GetPos_2() );
                        AddExtra( nq, n_q );
                    }
                    while( bui.Next() );
                }
            }
        }

        // this is a helper function for NextSubPass(); add 'pqn' to
        // 'marked_list' and increase extra by the size of the query
        // data structure plus raw query data
        //
        void SetupMateMark( 
                TQNum pqn, size_t & extra, 
                std::vector< TQNum > & marked_list )
        {
            if( !this->queries_.IsMarked( pqn ) ) {
                marked_list.push_back( pqn );
                size_t n_words( 3 + ((this->queries_.Len( pqn ) - 1)>>WSHIFT) );
                n_words <<= 1;
                extra += sizeof( CQueryData ) + sizeof( TWord )*n_words;
            }
        }

        // function object to set up query data structure for hash 2 pass
        //
        struct HashSetup
        {
            void operator()( CQueryData & qh ) const
            {
                TPrefix prefix( qh.Prefix() ), rc_prefix( 0 );
                ReverseComplement< SEQDATA_CODING >( rc_prefix, prefix );
                qh.SetPrefix( (prefix <= rc_prefix) ? prefix : rc_prefix );
                qh.SetStrand( (prefix < rc_prefix) ? STRAND_FW : STRAND_RV );
                qh.SetPalindrome( prefix == rc_prefix );
                qh.SetCurrHashIdx( 3 );
                qh.SetRightExtDir( qh.RightExtDir() );
                qh.SetExtLen();
            }
        };

        // add a single result of query expansion to the query data
        // cache
        //
        void AddExtra( const CQueryData & q, size_t & n_q )
        {
            qcache_.push_back( q );
            HashSetup()( *qcache_.rbegin() );
            ++n_q;
        }

        // instance constructor
        // 
        CSearchPass_ByHash( const CSearchPassDef::SInitData & init_data )
            : TBase( init_data ), start_idx_( 0 ), end_idx_( 0 ),
              wbuf_( 4*1024, 0 )
        { 
            typedef typename TBase::CException TException;

            if( MIN_ERR > this->n_err_ ) {
                M_THROW( TException, PASS_SKIP, "" ); 
            }

            qcache_.reserve( INIT_QCACHE_SIZE );
            NextSubPass();
        }

        // prepare for the next sub-pass (no additional sub-passes for hash 1)
        //
        void NextSubPass( void );

    private:

        //######################################################################
        //
        // SUBPASS HANDLING DATA
        //
        //######################################################################

        // helper constants: base letters in the query word and it binary log
        //
        static const size_t WLETTERS = 
            sizeof( TWord )*SCodingTraits< SEQDATA_CODING >::PACK_FACTOR;
        static const size_t WSHIFT = SBinLog< WLETTERS >::VALUE;

        // initial size of the temporary query data storage
        //
        static const size_t INIT_QCACHE_SIZE = 1024;

        // specializes NextSubPass() based on whether the search is paired
        //
        void NextSubPass_Ext( void );

        size_t start_idx_, // start and end indexes (relative to batch start)
               end_idx_;   // of queries searched in the current subpass

        std::vector< CQueryData > qcache_; // temporary storage of query data
                                           //   resulting from expansion of
                                           //   one query
        std::vector< TWord > wbuf_;        // storage for the raw data of the
                                           //   the query being expanded
};

//------------------------------------------------------------------------------
// CSearchPass functionality that must by specialized by 'paired'
//
template< int search_mode, int hash, bool paired > class CSearchPass_ByPaired;

//##############################################################################
template< int search_mode, int hash >
class CSearchPass_ByPaired< search_mode, hash, true > 
    : public virtual CSearchPass_Base< search_mode >
{
    private:

        typedef CSearchPass_Base< search_mode > TBase;

        // initial size of the subject segment list of in-place aligner
        static const size_t INIT_NUM_SEGS = 2;

    protected:

        typedef typename TBase::TPos TPos;

        // representation of a subject segment: the first member is the
        // start of the position in the database; the second member is
        // the length of the segment in bases
        //
        typedef std::pair< TPos, TSeqSize > TSegDescr;

        // query iterator type to use by in-place aligner
        //
        typedef CSeqFwIterator< SEQDATA_CODING, TWord, Uint1 > TQIter;

        // the type of in-place aligner
        //
        typedef CInPlaceAlign TAligner;

        // instance constructor
        CSearchPass_ByPaired( const CSearchPassDef::SInitData & init_data )
            : TBase( init_data ), ipam_vec_( init_data.ipam_vec )
        {
            results_.reserve( this->queries_.ResLimit() );
            ip_segs_.reserve( INIT_NUM_SEGS );
        }

        // check if in-place search is needed for a unique query
        //
        bool NeedsSearch( const CQueryData & q, TQNum qn ) const;

        // perform in-place search for queries in [qs, qe)
        //
        void SearchInPlace( 
                const CHit & hit, TQNum * qs, TQNum * qe, 
                TPos p, int n_err, TAligner & aligner );

        // adjust left and right segments for in-place matcher by
        // mate read length
        //
        void AdjustInPlaceSegs( 
                TSeqSize qlen, TPos anchor_start, TPos anchor_end ) {
            for( TAligner::TSegs::iterator i( ip_segs_.begin() );
                    i != ip_segs_.end(); ++i ) {
                if( i->left ) {
                    TSeqSize l( i->len + qlen );

                    if( i->pos + l > anchor_end ) i->len = anchor_end - i->pos;
                    else i->len = l;
                }
                else {
                    if( i->pos < qlen + anchor_start ) {
                        i->len += (i->pos - anchor_start);
                        i->pos = anchor_start;
                    }
                    else { i->pos -= qlen; i->len += qlen; }
                }
            }
        }

        // perform in-place alignment for the mates of the query and
        // its duplicates, and save the results; if no paired alignment
        // is found for a query pair, then the unpaired result is saved 
        // in some cases
        //
        void PostProcessMatch( const CHit & hit, int n_err );

        S_IPAM ipam_vec_;               // paired result configuration spec
        TAligner::TResults results_;    // temporary results storage for 
                                        //      in-place aligner
        TAligner::TSegs    ip_segs_;    // subject segment list for in-place
                                        //      aligner
};

//##############################################################################
template< int search_mode, int hash >
class CSearchPass_ByPaired< search_mode, hash, false > 
    : public virtual CSearchPass_Base< search_mode >
{
    private:

        typedef CSearchPass_Base< search_mode > TBase;

    protected:

        typedef typename TBase::TPos TPos;

        CSearchPass_ByPaired( const CSearchPassDef::SInitData & init_data )
            : TBase( init_data )
        {}

        // Save the result if necessary
        //
        void PostProcessMatch( const CHit & hit, int );
};

//------------------------------------------------------------------------------
// class performing one search pass (with possible multiple sub-passes if
// queries are being expanded)
//
template< int search_mode, int hash, bool paired >
class CSearchPass : public CSearchPass_ByHash< search_mode, hash, paired >,
                    public CSearchPass_ByPaired< search_mode, hash, paired >
{
    private:

        typedef CSearchPass_Base< search_mode > TBase;
        typedef CSearchPass_ByHash< search_mode, hash, paired > TBaseByHash;
        typedef CSearchPass_ByPaired< search_mode, hash, paired > TBaseByPaired;

        typedef typename TBase::TPos TPos;

    public:

        // instance constructor
        //
        CSearchPass( const CSearchPassDef::SInitData & init_data )
            : TBase( init_data ), 
              TBaseByHash( init_data ), TBaseByPaired( init_data ),
              skip_good_( false ), skip_bad_( false )
        {}

        // search the current batch of queries for the given hash value
        //
        void Run( void );

        // do not look at 'good' hashes in this pass
        //
        void SkipGood( void ) { skip_good_ = true; }

        // do not look at 'bad' hashes in this pass
        //
        void SkipBad( void ) { skip_bad_ = true; }

    private:

        // process a single pair of a query and a subject position that 
        // match on a 16-mer hash value
        //
        void ProcessCandidate( 
                const CQueryData & q, TPos pos, int n_err, TStrand s ); 

        // process one query (really a class of duplicate queries) 
        // corresponding to a "good" 16-mer hash value.
        //
        void ProcessQuery( const CQueryData & query );

        // process block of queries corresponding to the same "good" 16-mer
        // hash value
        //
        void ProcessQueryBlock( 
                const CQueryData & qstart, 
                const CQueryData & qend
                );

        // process block of special queries for a fixed subject position
        //
        void ProcessSpecialQueryBlockPos(
                CQueryData * qbase,
                typename TBase::TQExtDataTable & edt,
                const typename TBase::TBNF::TQExtSet & q_ext_set,
                TPos pos, int n_err, TStrand sstrand );

        // process block of queries corresponding to the same "bad"
        // 16-mer for a specific parameters of the number of errors,
        // query extension length, and query/subject strand configuration
        //
        void ProcessSpecialQueryBlock_Ext(
                CQueryData & qstart,
                CQueryData & qend,
                int n_err, int n_err_filter, TSeqSize ext_len,
                TStrand qstrand, TStrand estrand );
        
        // process block of unique query extensions against a
        // block of subject extensions.
        //
        void ProcessSpecialQueryBlock_ExtExt(
                CQueryData * qbase,
                const CIndexIterator::TExtData & ext,
                typename TBase::TQExtDataTable & edt,
                int n_err, int n_err_filter, TSeqSize ext_len,
                TStrand qstrand, TStrand estrand );

        // process block of queries corresponding to the same "bad" 16-mer
        // hash value and the same current hash index
        //
        void ProcessSpecialQueryBlock( 
                CQueryData & qstart, 
                CQueryData & qend );

        bool skip_good_;            // skip good n_mers
        bool skip_bad_;             // skip bad n_mers
};

END_NS( srprism )
END_STD_SCOPES

#include "search_pass_priv.hpp"

#endif

