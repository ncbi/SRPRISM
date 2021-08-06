/*  $Id$
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
 * File Description: seqinput template class going over multiple streams
 *                   serially
 *
 */

#ifndef __AM_SEQ_SERIAL_STREAM_HPP__
#define __AM_SEQ_SERIAL_STREAM_HPP__

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
struct CSerialStreamException : public common::CException
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

    M_EXCEPT_CTOR( CSerialStreamException )
};

//------------------------------------------------------------------------------
template< typename t_stream >
class CSerialStream : public CSeqInput
{
    private:

        static const int N_COLS = 1;

        typedef t_stream TStream;
        typedef TData::TSeq TSeq;

    public:

        CSerialStream( 
                const std::string & name, common::CFileBase::TCompression c )
            : d_( seq_, 0 ), c_( c )
        {
            done_ = true;

            if( !name.empty() )
            {
                for( std::string::size_type posb( 0 ), pose( posb );
                     pose != std::string::npos;
                     posb = pose + 1 )
                {
                    pose = name.find( ",", posb );
                    names_.push_back( name.substr( posb, pose - posb ) );
                }

                SRPRISM_ASSERT( !names_.empty() );
                stream_.reset( new TStream( names_[stream_idx_++], c_ ) );
                done_ = false;
            }
        }

        virtual int NCols( void ) const { return N_COLS; }

        virtual bool Next( void )
        {
            if( done_ ) return false;

            while( !stream_->Next() )
            {
                if( stream_idx_ >= names_.size() ) return false;
                stream_.reset( new TStream( names_[stream_idx_++], c_ ) );
            }

            seq_ = stream_->Data().seq;
            qual_ = stream_->Qual();
            id_ = stream_->Id() + '.' + std::to_string( stream_idx_ );
            title_ = stream_->Title();
            d_.size = stream_->Data().size;
            done_ = (stream_idx_ >= names_.size() && stream_->Done());
            return true;
        }

        virtual const TData & Data( int col ) const
        {
            if( col == 0 ) return d_;
            else M_THROW(
                    CSerialStreamException, COL, "invalid column: " << col );
        }

        virtual const TQual & Qual( int col ) const { return qual_; }

    private:

        std::vector< std::string > names_;
        std::unique_ptr< TStream > stream_;
        size_t stream_idx_ = 0;
        TSeq seq_;
        TQual qual_;
        TData d_;
        common::CFileBase::TCompression c_;
};

//------------------------------------------------------------------------------
inline std::auto_ptr< CSeqInput > MakeSerialStream( 
        const std::string & type, const std::string & name,
        common::CFileBase::TCompression c )
{
    CStreamFactory::TStreamType typeval( CStreamFactory::StreamType( type ) );

    switch( typeval ) {
        case CStreamFactory::STREAM_TYPE_FASTA: 
            return std::auto_ptr< CSeqInput >( 
                    new CSerialStream< CFastaStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_FASTQ:
            return std::auto_ptr< CSeqInput >( 
                    new CSerialStream< CFastqStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_CFASTA:
            return std::auto_ptr< CSeqInput >( 
                    new CSerialStream< CColorFastaStream >( name, c ) );
        case CStreamFactory::STREAM_TYPE_CFASTQ:
            return std::auto_ptr< CSeqInput >( 
                    new CSerialStream< CColorFastqStream >( name, c ) );
        default: M_THROW( CSerialStreamException, TYPE, "(" << name << ")" );
    }
}

END_NS( seq )
END_STD_SCOPES

#endif

