/*  $Id: seqinput_multistream.cpp 319730 2011-07-25 15:03:06Z morgulis $
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
 * File Description: input parser that reads multiple input streams
 *                   in parallel as columns
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <cctype>
#include <iostream>

#include "stream_factory.hpp"
#include "seqinput_multistream.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
CSeqInputMultiStream::CSeqInputMultiStream(
        const std::string & name, const std::string & type, int max_cols,
        CFileBase::TCompression c )
    : seq_n_( 0 )
{
    done_ = false;

    if( name.empty() ) {
        streams_.push_back( 
                CStreamFactory::MakeSeqStream( type, "", c ).release() );
    }
    else {
        for( std::string::size_type spos( 0 ), epos( 0 );
                epos != std::string::npos &&
                streams_.size() < (size_t)max_cols; ) {
            epos = name.find_first_of( ",", spos );
            std::string fname( name.substr( spos, epos - spos ) );
            streams_.push_back( CStreamFactory::MakeSeqStream( 
                        type, fname, c ).release() );
            spos = epos + 1;
        }
    }

    if( NCols() == 0 ) done_ = true;
}

//------------------------------------------------------------------------------
CSeqInputMultiStream::~CSeqInputMultiStream(void)
{
    typedef TStreams::iterator TIter;
    for( TIter i = streams_.begin(); i != streams_.end(); ++i ) delete *i;
}

//------------------------------------------------------------------------------
bool CSeqInputMultiStream::Next( void )
{
    if( done_ ) M_THROW( CException, EOS, "at " << seq_n_ );
    if( !streams_[0]->Next() ) { done_ = true; return false; }
    title_ = streams_[0]->Title();
    
    if( streams_.size() < 2 ) id_ = streams_[0]->Id();
    else {
        for( size_t j( 1 ); j < streams_.size(); ++j ) {
            if( !streams_[j]->Next() ) {
                M_THROW( CException, COL,
                         "premature end of stream at sequence " << seq_n_ <<
                         ", column " << j );
            }
        }

        const TSeqId & id_0( streams_[0]->Id() );
        size_t min_sz( id_0.size() );

        for( size_t j( 1 ); j < streams_.size(); ++j ) {
            min_sz = std::min( min_sz, streams_[j]->Id().size() );
        }

        size_t sz( min_sz );

        for( size_t j( 1 ); j < streams_.size(); ++j ) {
            const TSeqId & id_j( streams_[j]->Id() );
            size_t k( 0 );
            for( ; k < min_sz; ++k ) if( id_0[k] != id_j[k] ) break;
            sz = std::min( sz, k );
        }

        bool match( false );

        // check for special case "[./_]d+"
        //
        if( sz > 0 && ( 
                    id_0[sz-1] == '/' || 
                    id_0[sz-1] == '.' || 
                    id_0[sz-1] == '_' ) ) {
            match = true;

            for( size_t j( 0 ); j < streams_.size() && match; ++j ) {
                const TSeqId & id_j( streams_[j]->Id() );

                if( streams_[j]->Id().size() > sz ) {
                    for( size_t k( sz ); k < id_j.size(); ++k ) {
                        if( !isdigit( id_j[k] ) ) { match = false; break; }
                    }

                    size_t id_num( (size_t)atoi( id_j.substr( sz ).c_str() ) );
                    if( id_num != j + 1 ) { match = false; break; }
                }
                else match = false;
            }

            if( match ) --sz;
        }

        if( !match ) {
            for( size_t j( 0 ); j < streams_.size(); ++j ) {
                if( sz < streams_[j]->Id().size() ) {
                    M_THROW( CException, COL, "ids do not match" );
                }
            }
        }

        if( sz == 0 ) {
            M_THROW( CException, COL, 
                     "ids do not match at sequence " << seq_n_ );
        }

        id_ = id_0.substr( 0, sz );
    }

    done_= streams_[0]->Done();
    ++seq_n_;
    return true;
}

END_NS( seq )
END_STD_SCOPES

