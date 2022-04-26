/*  $Id: paired_stream.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: seqinput template class wrapping seqstreams containing
 *                   both mates in a single file
 *
 */

#ifndef __AM_SEQ_PAIRED_STREAM_HPP__
#define __AM_SEQ_PAIRED_STREAM_HPP__

#include "../common/def.h"

#include <string>

#ifndef NCBI_CPP_TK
#   include <seq/seqinput.hpp>
#   include <seq/stream_factory.hpp>
#   include <seq/fasta_stream.hpp>
#   include <seq/fastq_stream.hpp>
#   include <seq/cfasta_stream.hpp>
#   include <seq/cfastq_stream.hpp>
#else
#   include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/stream_factory.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/fasta_stream.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/fastq_stream.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/cfasta_stream.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/cfastq_stream.hpp>
#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
struct CPairedStreamException : public common::CException
{
    typedef common::CException TBase;

    static const TErrorCode COL  = 1;
    static const TErrorCode TYPE = 2;

    virtual const std::string ErrorMessage( TErrorCode code ) const
    {
        if( code == COL ) return "inconsistent columns";
        else if( code == TYPE ) return "bad stream type";
        else return TBase::ErrorMessage( code );
    }

    M_EXCEPT_CTOR( CPairedStreamException )
};

//------------------------------------------------------------------------------
template< typename t_stream >
class CPairedStream : public CSeqInput
{
    private:

        static const int N_COLS = 2;

        typedef t_stream TStream;
        typedef TData::TSeq TSeq;

    public:

        CPairedStream( 
                const std::string & name, common::CFileBase::TCompression c )
            : stream_( name, c ), 
              d0_( seq_[0], 0 ), d1_( seq_[1], 0 )
        {
            done_ = false;
        }

        virtual int NCols( void ) const { return N_COLS; }

        virtual bool Next( void )
        {
            if( !stream_.Next() ) return false;
            seq_[0] = stream_.Data().seq;
            qual_[0] = stream_.Qual();
            sz_[0] = stream_.Data().size;
            id_[0] = stream_.Id();
            title_[0] = stream_.Title();

            if( !stream_.Next() ) {
                M_THROW( CPairedStreamException, COL, "odd number of rows" );
            }

            seq_[1] = stream_.Data().seq;
            qual_[1] = stream_.Qual();
            sz_[1] = stream_.Data().size;
            id_[1] = stream_.Id();
            title_[1] = stream_.Title();
            size_t min_sz( std::min( id_[0].size(), id_[1].size() ) ),
                   sz( 0 );
            for( ; sz < min_sz; ++sz ) if( id_[0][sz] != id_[1][sz] ) break;
            bool match( false );

            if( sz > 0 && (
                    id_[0][sz-1] == '/' ||
                    id_[0][sz-1] == '.' ||
                    id_[0][sz-1] == '_' ) ) {
                match = true;

                for( size_t j( 0 ); j < (size_t)N_COLS; ++j ) {
                    if( match && id_[j].size() > sz ) {
                        for( size_t k( sz ); k < id_[j].size(); ++k ) {
                            if( !isdigit( id_[j][k] ) ) { 
                                match = false; break; 
                            }
                        }

                        if( match ) {
                            size_t id_num( (size_t)atoi( 
                                        id_[j].substr( sz ).c_str() ) );
                            if( id_num != j + 1 ) { match = false; break; }
                        }
                    }
                    else match = false;
                }

                if( match ) --sz;
            }

            if( !match ) {
                for( size_t j( 0 ); j < (size_t)N_COLS; ++j ) {
                    if( sz < id_[j].size() ) {
                        M_THROW( CPairedStreamException, COL, 
                                 "ids do not match: " <<
                                 id_[0] << " , " << id_[1] );
                    }
                }
            }

            if( sz == 0 ) {
                M_THROW( CPairedStreamException, COL, 
                         "ids do not match: " << id_[0] << " , " << id_[1] );
            }

            CSeqInput::id_ = id_[0].substr( 0, sz );
            done_ = stream_.Done();
            d0_.size = sz_[0];
            d1_.size = sz_[1];
            return true;
        }

        virtual const TData & Data( int col ) const 
        { 
            if( col == 0 ) return d0_;
            else if( col == 1 ) return d1_;
            else M_THROW( 
                    CPairedStreamException, COL, "invalud column: " << col );
        }

        virtual const TQual & Qual( int col ) const { return qual_[col]; }

    private:

        TStream stream_;
        TSeq seq_[N_COLS];
        TQual qual_[N_COLS];
        TSeqSize sz_[N_COLS];
        TSeqId id_[N_COLS];
        TSeqTitle title_[N_COLS];
        TData d0_, d1_;
};

//------------------------------------------------------------------------------
// inline std::auto_ptr< CSeqInput > MakePairedStream( 
inline std::unique_ptr< CSeqInput > MakePairedStream( 
        const std::string & type, const std::string & name,
        common::CFileBase::TCompression c )
{
    CStreamFactory::TStreamType typeval( CStreamFactory::StreamType( type ) );

    switch( typeval ) {
        case CStreamFactory::STREAM_TYPE_FASTA: 
            // return std::auto_ptr< CSeqInput >( 
            return std::unique_ptr< CSeqInput >( 
                    new CPairedStream< CFastaStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_FASTQ:
            // return std::auto_ptr< CSeqInput >( 
            return std::unique_ptr< CSeqInput >( 
                    new CPairedStream< CFastqStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_CFASTA:
            // return std::auto_ptr< CSeqInput >( 
            return std::unique_ptr< CSeqInput >( 
                    new CPairedStream< CColorFastaStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_CFASTQ:
            // return std::auto_ptr< CSeqInput >( 
            return std::unique_ptr< CSeqInput >( 
                    new CPairedStream< CColorFastqStream >( name, c ) );
        default: M_THROW( CPairedStreamException, TYPE, "(" << name << ")" );
    }
}

END_NS( seq )
END_STD_SCOPES

#endif

