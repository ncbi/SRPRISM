/*  $Id: fastq_stream.cpp 342358 2011-10-27 14:24:39Z morgulis $
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
 * File Description: fastq parser implementation
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <cassert>

#include "fastq_stream.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
CFastqStream::CFastqStream( 
        const std::string & name, CFileBase::TCompression c )
    : CStreamBase( name, c ), state_( E_START )
{ ReadLine(); }

//------------------------------------------------------------------------------
bool CFastqStream::ReadLine(void)
{
    if( in_->Eof() ) { done_ = true; return false; }
    else {
        try { line_ = in_->GetLine(); }
        catch( CReadTextFile::CException & ) {}

        return true;
    }
}

//------------------------------------------------------------------------------
void CFastqStream::ProcessDefline(void)
{
    const char * SEPS = " \t";
    std::string::size_type pos = line_.find_first_of( SEPS, 1 );
    id_ = line_.substr( 1, pos - 1 );
    
    if( pos != std::string::npos ) {
        pos = line_.find_first_not_of( SEPS, ++pos );
        if( pos != std::string::npos ) title_ = line_.substr( pos );
    }
}

//------------------------------------------------------------------------------
void CFastqStream::ProcessQTitle( void )
{
    const char * SEPS = " \t";
    std::string::size_type pos = line_.find_first_of( SEPS, 1 );
    TSeqId id = line_.substr( 1, pos - 1 );

    if( !id.empty() && id != id_ ) {
        M_THROW( CException, PARSE,
                 "id mismatch at "
                 << name_ << ":" << in_->LineNo() );
    }
    
    if( pos != std::string::npos ) {
        pos = line_.find_first_not_of( SEPS, ++pos );
        TSeqTitle title = line_.substr( pos );
        
        if( title != title_ ) {
            M_THROW( CException, PARSE,
                     "title mismatch at "
                     << name_ << ":" << in_->LineNo() );
        }
    }
}

//------------------------------------------------------------------------------
void CFastqStream::ProcessLine(void)
{
    switch( state_ ) {
        case E_START:
            seq_.clear(); qual_.clear(); id_.clear(); title_.clear();

            if( line_.empty() || line_[0] == '#' ) { 
                if( ! ReadLine() ) state_ = E_EMPTY;
                return;
            }

            if( line_[0] == '@' ) { ProcessDefline(); state_ = E_DATA; }
            else {
                M_THROW( CException, PARSE,
                         "expected '@' at " << 
                         name_ << ":" << in_->LineNo() );
            }

            ReadLine();
            break;

        case E_DATA:
            if( line_.empty() || line_[0] == '#' ) { ReadLine(); return; }

            if( line_[0] == '+' ) { 
                ProcessQTitle(); 
                ReadLine(); 
                state_ = E_QDATA;
            }
            else {
                if( line_.find_first_not_of( 
                            SCodingTraits< CODING >::ALPHABET_STRING ) !=
                        std::string::npos ) {
                    M_THROW( CException, PARSE,
                             "illegal letter at " << 
                             name_ << ":" << in_->LineNo() );
                }

                std::copy( line_.begin(), line_.end(), 
                           std::back_inserter( seq_ ) );
                ReadLine();
            }

            break;

        case E_QDATA:
            if( line_.empty() ) { ReadLine(); return; }
            std::copy( line_.begin(), line_.end(), std::back_inserter( qual_ ) );

            if( qual_.size() > seq_.size() ) {
                M_THROW( CException, PARSE,
                        "sizes of sequence data and quality data "
                        << "do not match at "
                        << name_ << ":" << in_->LineNo() );
            }

            if( qual_.size() >= seq_.size() ) {
                seq_qual_.clear();
                std::copy( 
                        qual_.begin(), qual_.end(), 
                        std::back_inserter( seq_qual_ ) );
                state_ = E_READY;

                do{ ReadLine(); }
                while( !done_ && (line_.empty() || line_[0] == '#') );
            }
            else ReadLine();

            break;

        default: SRPRISM_ASSERT( false );
    }
}

//------------------------------------------------------------------------------
bool CFastqStream::Next(void)
{
    if( done_ ) return false;
    while( !done_ && state_ != E_READY && state_ != E_EMPTY ) ProcessLine();
    if( state_ == E_EMPTY ) return false;

    if( state_ == E_QDATA ) {
        if( seq_.size() != qual_.size() ) {
            M_THROW( CException, PARSE,
                     "sizes of sequence data and quality data "
                     << "do not match at "
                     << name_ << ":" << in_->LineNo() );
        }

        state_ = E_READY;
    }

    if( state_ != E_READY ) {
        M_THROW( CException, PARSE,
                 "unexpected end of file at " 
                 << name_ << ":" << in_->LineNo() );
    }

    seq_data_.size = seq_.size();
    state_ = E_START;
    return true;
}

END_NS( seq )
END_STD_SCOPES
