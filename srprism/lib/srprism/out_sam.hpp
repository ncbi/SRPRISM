/*  $Id: out_sam.hpp 536751 2017-05-23 13:07:55Z morgulis $
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
 * File Description: output in SAM format
 *
 */

#ifndef __SRPRISM_OUT_SAM_HPP__
#define __SRPRISM_OUT_SAM_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <common/exception.hpp>

#include <srprism/out_base.hpp>
#include <srprism/result.hpp>

#else

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>

#include <../src/internal/align_toolbox/srprism/lib/srprism/out_base.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class COutSAM : public COutBase
{
    private:

        static const std::string SEP;

        struct SSAMRecord
        {
            typedef unsigned long TFlags;

            static const TFlags PAIRED_QUERY_FLAG  = 0x1;
            static const TFlags PAIRED_ALIGN_FLAG  = 0x2;
            static const TFlags SEQ_UNMAPPED_FLAG  = 0x4;
            static const TFlags MATE_UNMAPPED_FLAG = 0x8;
            static const TFlags SEQ_STRAND_FLAG    = 0x10;
            static const TFlags MATE_STRAND_FLAG   = 0x20;
            static const TFlags SEQ_FIRST_FLAG     = 0x40;
            static const TFlags SEQ_SECOND_FLAG    = 0x80;
            static const TFlags NOT_PRIMARY_FLAG   = 0x100;

            struct STag
            {
                static const int TAG_NAME_LEN = 2;

                static const char * VALID_TYPES;

                char name[TAG_NAME_LEN];
                char type;
                std::string value;

                STag( const char * n, char type, const std::string & value )
                    : type( type ), value( value )
                { for( int i = 0; i < TAG_NAME_LEN; ++i ) name[i] = n[i]; }
            };

            typedef std::vector< STag > TTags;

            std::string qname;
            TFlags flags;
            std::string sname;
            std::string msname;
            common::Sint8 pos;
            std::string CIGAR_str;
            common::Sint8 mpos;
            common::Sint8 diff;
            std::string seq;
            TTags tags;
            int quality;
            std::string qstr;
            std::string extra_tags;

            SSAMRecord( void )
                : flags( 0 ), pos( 0 ), mpos( 0 )
            {}

            void AddTag( const char * name, char type, const std::string & value )
            { tags.push_back( STag( name, type, value ) ); }

            void AddITag( const char * name, char type, common::Sint8 value )
            {
                std::ostringstream os;
                os << value;
                AddTag( name, type, os.str() );
            }

            const std::string Format( void ) const
            {
                static const TFlags pflags( 
                        PAIRED_QUERY_FLAG|PAIRED_ALIGN_FLAG );
                bool paired( (flags&pflags) == pflags );

                common::Sint8 d( diff );

                if( paired ) {
                    if( mpos - pos < 0 ) d = -diff;
                    else if( mpos == pos ) {
                        if( (flags&0x10) == 0 ) {
                            if( (flags&0x20) == 0 ) {
                                if( (flags&0x40) == 0 ) d = -diff;
                            }
                        }
                        else {
                            if( (flags&0x20) == 0 ) d = -diff;
                            else if( (flags&0x40) == 0 ) d = -diff;
                        }
                    }
                }

                std::ostringstream os;
                os << qname      << '\t'
                   << flags      << '\t'
                   << sname      << '\t'
                   << pos        << '\t'
                   << quality    << '\t'
                   << CIGAR_str  << '\t'
                   << (msname.empty() ? std::string( "*" ) : 
                        (msname == sname ? std::string( "=" ) : msname)) << '\t'
                   << mpos       << '\t'
                   << (paired ? d : 0) << '\t'
                   << seq << "\t" << qstr;
            
                for( TTags::const_iterator i( tags.begin() ); 
                        i != tags.end(); ++i ) {
                    os << '\t' << i->name[0] << i->name[1] << ':'
                       << i->type << ':' << i->value;
                }

                if( !extra_tags.empty() ) {
                    os << '\t' << extra_tags;
                }

                return os.str();
            }
        };

    public:

        COutSAM( 
                const std::string & name,
                const std::string & input,
                const std::string & input_fmt,
                const std::string & extra_tags,
                const std::string & cmdline,
                bool print_header,
                common::CFileBase::TCompression input_c, 
                bool skip_unmapped,
                bool force_paired,
                bool force_unpaired,
                bool no_qids,
                bool out_xa,
                CSeqStore * seq_store,
                CSIdMap * sid_map )
            : COutBase( 
                    name, input, input_fmt, input_c, 
                    skip_unmapped, force_paired, force_unpaired, no_qids,
                    seq_store, sid_map ),
              extra_tags_( extra_tags ),
              out_xa_( out_xa )
        {
            if( print_header )
            {
                (*os_) << "@HD\tVN:1.0\tGO:query\n";
                (*os_) << "@PG\tID:srprism\tPN:srprism\tCL:" << cmdline << '\n';
                seq_store->Load();

                for( size_t i( 0 ); i < seq_store->NSeq(); ++i )
                {
                    (*os_) << "@SQ\tSN:" << (*sid_map)[i]
                           << "\tLN:" << seq_store_->GetSeqLen( i ) << '\n';
                }

                seq_store->Unload();
            }
        }

        ~COutSAM();

        void EmptyOut( 
                int idx, common::Uint4 mpos = 0, 
                const std::string & sname = "",
                seq::TStrand mstrand = seq::STRAND_FW,
                bool primary = false );

        virtual void FinalizeBatch();

    private:

        static std::string ComputeCIGAR( const CResult & result, int idx );
        std::string ComputeMDTag( const CResult & result, int idx );
        std::string ComputeSeq( int idx, bool reverse );
        void SetUpQData( int idx, SSAMRecord & r );

        virtual void ResultOut( 
                const CResult & result, bool mate_unmapped, 
                TQueryOrdId q_adj, int const * pg, bool primary = true );

        virtual void ResultOut(
                const CResult & result_1,
                const CResult & result_2,
                TQueryOrdId q_adj, int const * pg );

        COutSAM( const COutSAM & );
        COutSAM & operator=( const COutSAM & );

        std::string extra_tags_;
        bool out_xa_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

