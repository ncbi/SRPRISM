/*  $Id: search.hpp 536751 2017-05-23 13:07:55Z morgulis $
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
 * File Description: representation of a search task
 *
 */

#ifndef __SRPRIZM_SEARCH_HPP__
#define __SRPRIZM_SEARCH_HPP__

#include "../common/def.h"

#include <string>
#include <memory>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/stat.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/seqstore.hpp>
#include <srprism/sidmap.hpp>
#include <srprism/batch.hpp>
#include <srprism/out_base.hpp>
#include <srprism/out_sam.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/stat.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/sidmap.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/batch.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/out_base.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CSearch
{
    public:

        static const common::Uint4 MIN_RES_LIMIT = 1UL;
        static const common::Uint4 MAX_RES_LIMIT = 254UL;
        static const common::Uint2 MAX_PAIR_FUZZ = 10*common::KILOBYTE;
        static const common::Uint2 MIN_QLEN      = MIN_QUERY_LEN;
        static const common::Uint2 MAX_QLEN      = MAX_QUERY_LEN;
        static const common::Uint1 MAX_N_ERR     = MAX_ERR;

        struct SOptions
        {
            SOptions()
                : input_fmt( "fasta" ),
                  tmpdir( "." ),
                  resconf_str( "0100" ),
                  input_compression( common::CFileBase::COMPRESSION_AUTO ),
                  mem_limit( 2048 ),
                  batch_limit( 10000000UL ),
                  start_batch( 1 ), end_batch( 1 ),
                  res_limit( 10 ),
                  repeat_threshold( 4096 ),
                  pair_distance( 500 ),
                  pair_fuzz( 490 ),
                  max_qlen( 16 ),
                  n_threads( 1 ),
                  sa_start( 1 ),
                  sa_end( 8192 ),
                  n_err( 0 ),
                  search_mode( srprism::SSearchMode::DEFAULT ),
                  force_paired( false ),
                  force_unpaired( true ),
                  use_sids( true ),
                  use_qids( true ),
                  skip_unmapped( true ),
                  strict_batch( false ),
                  discover_sep( false ),
                  discover_sep_stop( false ),
                  randomize( false ),
                  random_seed( false ),
                  use_fixed_hc( false )
            {
            }

            std::string index_basename;
            std::string input;
            std::string output;
            std::string input_fmt;
            std::string tmpdir;
            std::string resconf_str;
            std::string paired_log;
            std::string extra_tags;
            std::string hist_fname;
            std::string cmdline;
            common::CFileBase::TCompression input_compression;
            size_t mem_limit;
            common::Uint8 batch_limit;
            common::Uint4 start_batch;
            common::Uint4 end_batch;
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
            bool force_paired;
            bool force_unpaired;
            bool use_sids;
            bool use_qids;
            bool skip_unmapped;
            bool strict_batch;
            bool discover_sep;
            bool discover_sep_stop;
            bool randomize;
            bool random_seed;
            bool use_fixed_hc;
            bool sam_header;
        };

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode VALIDATE = 0;
            static const TErrorCode INPUT = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == VALIDATE ) return "validation error";
                else if( code == INPUT ) return "input error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CSearch( const SOptions & options );
        ~CSearch(void);
        void Run(void);

    private:

        CSearch( const CSearch & );
        CSearch & operator=( const CSearch & );

        void Validate( const SOptions & options ) const;
        void Run_priv(void);

        std::shared_ptr< CMemoryManager > mem_mgr_p_;
        std::auto_ptr< CSIdMap > sidmap_p_;
        std::auto_ptr< CSeqStore > seqstore_p_;

        std::unique_ptr< common::CTmpStore > tmp_store_p_;
        std::unique_ptr< COutSAM_Collator > out_p_;

        std::string input_;
        std::string input_fmt_;
        std::string extra_tags_;

        common::CFileBase::TCompression input_c_;

        bool use_sids_;
        bool force_paired_;
        bool force_unpaired_;
        bool strict_batch_;
        bool skip_unmapped_;
        bool use_qids_;

        Uint4 start_batch_;
        Uint4 end_batch_;
        Uint8 batch_limit_;

        CBatch::SBatchInitData batch_init_data_;
        CStatMap global_stats_;
};

END_NS( srprism )
END_STD_SCOPES


#endif

