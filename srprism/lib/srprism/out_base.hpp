/*  $Id: out_base.hpp 431273 2014-04-02 17:10:44Z morgulis $
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
 * File Description: base class for results output
 *
 */

#ifndef __SRPRISM_OUT_BASE_HPP__
#define __SRPRISM_OUT_BASE_HPP__

#include "../common/def.h"

#include <vector>

#ifndef NCBI_CPP_TK

#include <seq/seqinput_factory.hpp>
#include <seq/seqinput.hpp>
#include <srprism/result.hpp>
#include <srprism/sidmap.hpp>
#include <srprism/query_store.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/seq/seqinput_factory.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/sidmap.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_store.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class COutBase
{
    protected:

        class CInputAdapter
        {
            public:
    
                CInputAdapter() : qn_( 0 ) {}

                virtual ~CInputAdapter() {}

                virtual std::string Qual( size_t idx ) const = 0;
                virtual bool Skip( TQueryOrdId qn ) = 0;
                virtual int NCols( void ) const = 0;
                virtual seq::TSeqId Id( void ) const = 0;
                virtual std::string QData( size_t idx ) const = 0;
                virtual bool Done( void ) const = 0;

                TQueryOrdId QId( void ) const { return qn_; }

            protected:

                TQueryOrdId qn_;
        };

        class CInputAdapterSpec : public CInputAdapter
        {
            public:
    
                CInputAdapterSpec( std::auto_ptr< seq::CSeqInput > in ) 
                    : in_( in ) {}

                virtual std::string Qual( size_t idx ) const 
                { return in_->Qual( idx ); }

                virtual bool Skip( TQueryOrdId qn )
                { 
                    while( qn_ < qn ) { 
                        if( !in_->Next() ) return false;
                        ++qn_; 
                    } 

                    return true;
                }

                virtual int NCols( void ) const { return in_->NCols(); }
                virtual seq::TSeqId Id( void ) const { return in_->Id(); }

                virtual std::string QData( size_t idx ) const
                {
                    const std::vector< common::Uint1 > & seq( 
                            in_->Data( idx ).seq );
                    std::string result;
                    std::copy( 
                            seq.begin(), seq.end(), 
                            std::back_inserter( result ) );
                    return result;
                }

                virtual bool Done( void ) const { return in_->Done(); }

            private:

                CInputAdapterSpec( const CInputAdapterSpec & );
                CInputAdapterSpec & operator=( const CInputAdapterSpec & );
    
                std::auto_ptr< seq::CSeqInput > in_;
        };

    public:

        typedef std::vector< CResult * > TResults;

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode OPEN = 0;
            static const TErrorCode INFMT = 1;
            static const TErrorCode READ = 2;
            static const TErrorCode NOT_IMPLEMENTED = 3;

            const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == OPEN ) return "open error";
                else if( code == INFMT ) return "bad input format";
                else if( code == READ ) return "read error";
                else if( code == NOT_IMPLEMENTED ) {
                    return "operation is not implemented";
                }
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException );
        };

        COutBase( 
                const std::string & name,
                const std::string & input,
                const std::string & input_fmt,
                common::CFileBase::TCompression input_c,
                bool skip_unmapped,
                bool force_paired,
                bool force_unpaired,
                bool no_qids,
                CSeqStore * seq_store,
                CSIdMap * sid_map )
            : os_( 0 ), os_p_( 0 ), in_p_( 0 ), 
              seq_store_( seq_store ), sid_map_( sid_map ), 
              input_fmt_( input_fmt ),
              skip_unmapped_( skip_unmapped ), paired_( false ),
              no_qids_( no_qids )
        {
            if( name.empty() ) os_ = &std::cout;
            else {
                os_ = new std::ofstream( name.c_str() );

                if( !os_ || !os_->good() ) {
                    M_THROW( CException, OPEN, "[" << name << "]" );
                }

                os_p_.reset( os_ );
            }

            int n_cols( 0 );

            if( force_unpaired ) n_cols = 1;
            else if( force_paired ) n_cols = 2;

            in_p_.reset( new CInputAdapterSpec( 
                        seq::CSeqInputFactory::MakeSeqInput( 
                            input_fmt, input, n_cols, input_c ) ) );
            if( in_p_->NCols() > 1 ) paired_ = true;
            if( force_paired ) paired_ = true;
            if( force_unpaired ) paired_ = false;
        }

        virtual ~COutBase() {}

        void ResultsOut( const TResults & results, TQueryOrdId q_adj, 
                         int * const pg = 0 )
        {
            if( results.empty() ) return;
            size_t n_res[3] = {0, 0, 0};

            typedef TResults::const_iterator TIter;

            for( TIter i( results.begin() ); i < results.end(); ++i ) {
                SRPRISM_ASSERT( (*i)->QOrdId( q_adj, paired_ ) == 
                        results[0]->QOrdId( q_adj, paired_ ) );
                SRPRISM_ASSERT( (*i)->PairPos() < 2 );
                if( (*i)->Paired() ) ++n_res[0];
                else ++n_res[(*i)->PairPos() + 1];
            }

            /*
            This was the original interpretation of mat_unmapped flag:
            changed to reflect the state of individual results.

            bool mate_unmapped( 
                    n_res[0] == 0 && (n_res[1] == 0 || n_res[2] == 0) );
            */

            bool mate_unmapped( n_res[0] == 0 );

            if( n_res[1] != 0 && n_res[2] != 0 ) {
                TIter i( results.begin() ), j( i );
                int ippos( (*i)->PairPos() );
                for( ; j != results.end() && (*j)->PairPos() == ippos ; ++j );
                SRPRISM_ASSERT( j != results.end() );
                
                try { ResultOut( **i, **j, q_adj, pg ); }
                catch( CException & e ) {
                    if( e.ErrorCode() == CException::NOT_IMPLEMENTED ) {
                        ResultOut( **i, mate_unmapped, q_adj, pg, true );
                        ++i;

                        for( ; i != results.end(); ++i ) {
                            ResultOut( **i, mate_unmapped, q_adj, pg, false );
                        }

                        return;
                    }

                    throw;
                }

                ++i;

                for( ; i != results.end(); ++i ) {
                    if( i != j ) {
                        ResultOut( **i, mate_unmapped, q_adj, pg, false );
                    }
                }
            }
            else {
                TIter i( results.begin() );
                ResultOut( **i, mate_unmapped, q_adj, pg, true );
                ++i;

                for( ; i != results.end(); ++i ) {
                    ResultOut( **i, mate_unmapped, q_adj, pg, false );
                }
            }
        }

        void SetUpQueryInfo( CQueryStore const * qs, TQueryOrdId q_adj ) {
            qs_ = qs;
            q_adj_ = q_adj;
        }

        virtual void FinalizeBatch() {}

    protected:

        virtual void ResultOut( 
                const CResult & result, 
                bool mate_unmapped,
                TQueryOrdId q_adj,
                int const * pg,
                bool primary = true ) = 0;

        virtual void ResultOut(
                const CResult & result_1,
                const CResult & result_2,
                TQueryOrdId q_adj,
                int const * pg )
        {
            M_THROW( CException, NOT_IMPLEMENTED, "" );
        }

        std::ostream * os_;
        std::auto_ptr< std::ostream > os_p_;
        std::auto_ptr< CInputAdapter > in_p_;
        CSeqStore * seq_store_;
        CSIdMap * sid_map_;
        std::string input_fmt_;
        bool skip_unmapped_;
        bool paired_;
        bool no_qids_;
        CQueryStore const * qs_;
        TQueryOrdId q_adj_;

    private:

        COutBase( const COutBase & );
        COutBase & operator=( const COutBase & );
};

END_NS( srprism )
END_STD_SCOPES

#endif

