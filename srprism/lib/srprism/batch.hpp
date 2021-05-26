/*  $Id: batch.hpp 591182 2019-08-12 16:55:27Z morgulis $
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

#ifndef __SRPRISM_BATCH_HPP__
#define __SRPRISM_BATCH_HPP__

#include "../common/def.h"

#include <atomic>
#include <string>
#include <memory>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <common/trace.hpp>
#include <common/tmpstore.hpp>
#include <seq/seqinput.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/seqstore.hpp>
#include <srprism/tmpres_mgr.hpp>
#include <srprism/query_store.hpp>
#include <srprism/search_pass.hpp>
#include <srprism/out_base.hpp>
#include <srprism/rmap.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/trace.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/tmpstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/tmpres_mgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_store.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/search_pass.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/out_base.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/rmap.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CBatch
{
    public:

        static const size_t TMP_RES_BUF_SIZE = 1024*1024ULL;

        struct SBatchInitData
        {
            std::string index_basename;
            std::string tmpdir;
            std::string paired_log;
            std::string hist_fname;
            std::string resconf_str;
            common::Uint8 batch_limit;
            common::Uint4 res_limit;
            common::Uint4 repeat_threshold;
            common::Uint4 fixed_hc;
            common::Uint2 pair_distance;
            common::Uint2 pair_fuzz;
            common::Uint2 max_qlen;
            common::Uint2 n_threads;
            common::Sint2 sa_start;
            common::Sint2 sa_end;
            common::Uint1 n_err;
            int search_mode;
            bool use_qids;
            bool use_sids;
            bool paired;
            bool use_fixed_hc;
            bool discover_sep;
            bool discover_sep_stop;
            bool randomize;
            bool random_seed;

            S_IPAM ipam_vec;

            std::shared_ptr< CMemoryManager > mem_mgr_p;
            CSeqStore * seqstore_p;

            void * u_tmp_res_buf;
            size_t u_tmp_res_buf_size;

            void * p_tmp_res_buf;
            size_t p_tmp_res_buf_size;

            CStatMap * search_stats;
        };

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode POSTPROCESS = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == POSTPROCESS ) return "post process error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CBatch( SBatchInitData & init_data, 
                CSeqInput & in, TQueryOrdId start_qid, Uint4 batch_oid );
        
        ~CBatch() { if( queries_p_.get() != 0 ) queries_p_->FreeQueryData(); }

        template< int hash, bool paired > void 
        RunPass( bool skip_good, bool skip_bad );

        template< bool paired > bool Run(void);

        static void RunBatchPaired( CBatch * batch );
        static void RunBatchSingle( CBatch * batch );

        template< int search_mode > bool DiscoverInsertSize( void );
        template< int search_mode > bool InterProcess( void );
        template< int search_mode, bool paired > void PostProcess(void);

        template< int search_mode, bool paired > 
        CResult * RemoveDuplicates( CResult * s, CResult * e );

        template< int search_mode, bool paired > 
        CResult * RemoveDuplicatesForSubject( 
                CResult * s, CResult * e, CResult * o );

        TQueryOrdId StartQId(void) const { return start_qid_; }
        TQueryOrdId EndQId(void) const { return end_qid_; }

        void SetBatchOutput( COutBase * out_p )
        { out_p_.reset( out_p ); }

        std::string GetTmpName( std::string const & pfx )
        { return tmp_store_.Register( pfx ); }

    private:

        static const char * TMP_RES_FNAME;

        static const size_t MIN_RES_LIMIT = 10;

        CBatch( const CBatch & );
        CBatch & operator=( const CBatch & );

        template< int search_mode, bool paired >
        std::pair< TQNum, TQNum > ComputeQNumBounds( 
                TQNum start, size_t max_res ) const;

        typedef std::pair< CResult *, CResult * > TPairCandidate;
        typedef std::vector< TPairCandidate > TPairCandidates;

        //----------------------------------------------------------------------
        //  Types, values, and methods used for insert size discovery.
        //
        static const TSeqSize ISD_MIN_RANGE = 200;
        static const TSeqSize ISD_MAX_RANGE = 32*1024;
        static const TSeqSize ISD_MAX_INSERT = 1000000UL;
        static const Uint4 ISD_RES_LIMIT = 255;

        static const double ISD_ACCEPT_PER_SC_RATIO;
        static const double ISD_ACCEPT_RATIO;

        typedef std::map< Uint4, Uint8 > THistogram;

        void UpdateHistogramForSubject( 
                CResult * ls, CResult * le, CResult * rs, CResult * re, 
                THistogram * h, bool * scv );

        void UpdateHistogram( CResult * s, CResult * e, THistogram * h );
        bool ProcessHistogram( THistogram * h );
        void DumpHistogram( THistogram * h );

        bool ProcessSample( 
                TSeqSize & min, TSeqSize range, Uint8 total, 
                const THistogram & h );
        //----------------------------------------------------------------------

        template< int search_mode >
        bool CombineUnpairedResults( 
                CResult * s, CResult * e, 
                size_t & n_found, size_t & n_left, size_t & n_cancelled );

        template< int search_mode >
        void RetainUnpairedResults( CResult * s, CResult * e );

        template< int search_mode >
        void UpdateCounts( CResult * s, CResult * e );

        void MarkDone( CResult * s, CResult * e );
        void MarkUPRes( CResult * s, CResult * e );
        void ClearDone4Search( CResult * s, CResult * e );

        template< int search_mode > 
        size_t IdentifyPairs( CResult * s, size_t n_left, size_t n_right );

        int GetStrandConfig( CResult * l, CResult * r );
        bool CheckStrandConfig( CResult * l, CResult * r );

        size_t IdentifyPairsForSubject( 
                CResult * ls, CResult * le, 
                CResult * rs, CResult * re, TPairCandidates & pc );

        void CheckPair( 
                CResult * l, CResult * r, 
                TSeqSize e, size_t & res, TPairCandidates & pc );

        void DiversifySubjects( CResult * s, CResult * e );

        int ProvidesGuarantee( TQNum qn, int nerr ) {
            CQueryStore & qs( *queries_p_.get() );

            if( search_mode_ == SSearchMode::DEFAULT ||
                search_mode_ == SSearchMode::SUM_ERR ) {
                if( !qs.HasRepHashes( qn ) ) return 1;
                if( qs.HasRepBUHash( qn ) ) return 0;

                if( qs.IsLong( qn ) ) {
                    if( nerr < qs.GetComplHashCount( qn ) ) return 1;
                    return 0;
                }
                else if ( nerr == 0 ) {
                    return qs.GetComplHashCount( qn ) > 0 ? 1 : 0;
                }
                else return nerr + 1 < qs.GetComplHashCount( qn ) ? 1 : 0;
            }
            else return qs.HasRepHashes( qn ) ? 0 : 1;
        }

        SBatchInitData init_data_;
        common::CTmpStore tmp_store_;
        std::unique_ptr< CTmpResMgr > u_tmpres_mgr_;
        std::unique_ptr< CTmpResMgr > p_tmpres_mgr_;
        CTmpResMgr * tmpres_mgr_p_;
        CSeqStore & seqstore_;
        CRMap rmap_;
        bool use_sids_, use_qids_;
        int search_mode_;
        Uint4 res_limit_;
        Uint4 final_res_limit_;
        Uint4 batch_oid_;
        TQueryOrdId start_qid_;
        TQueryOrdId end_qid_;
        std::auto_ptr< CQueryStore > queries_p_;
        CSearchPassDef::SInitData pass_init_data_;
        std::unique_ptr< COutBase > out_p_;
        std::string paired_log_;

public:

        std::atomic< bool > done_;
};

END_NS( srprism )
END_STD_SCOPES

#include "batch_priv.hpp"

#endif

