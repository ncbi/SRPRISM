/*  $Id: fasta_stream.cpp 205515 2010-09-20 15:51:04Z morgulis $
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
 * File Description: fasta parser implementation
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include <cassert>

#include "fasta_stream.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
CFastaStream::CFastaStream( 
        const std::string & name, CFileBase::TCompression c )
    : CStreamBase( name, c ), state_( E_START )
{ ReadLine(); }

//------------------------------------------------------------------------------
bool CFastaStream::ReadLine(void)
{
    if( in_->Eof() ) { done_ = true; state_ = E_READY; return false; }
    else {
        try { line_.clear(); line_ = in_->GetLine(); }
        catch( CReadTextFile::CException & ) {}

        return true;
    }
}

//------------------------------------------------------------------------------
void CFastaStream::ProcessDefline(void)
{
    const char * SEPS = " \t";
    std::string::size_type pos = line_.find_first_of( SEPS, 1 );
    id_ = line_.substr( 1, pos - 1);
    
    if( pos != std::string::npos ) {
        pos = line_.find_first_not_of( SEPS, ++pos );
        if( pos != std::string::npos ) title_ = line_.substr( pos );
    }
}

//------------------------------------------------------------------------------
void CFastaStream::ProcessLine(void)
{
    switch( state_ ) {
        case E_START:
            seq_.clear(); id_.clear(); title_.clear();

            if( line_.empty() || line_[0] == '#' ) { 
                if( !ReadLine() ) state_ = E_EMPTY;
                return; 
            }

            if( line_[0] == '>' ) { ProcessDefline(); state_ = E_DATA; }
            else {
                M_THROW( CException, PARSE,
                         "expected defline at " << 
                         name_ << ":" << in_->LineNo() );
            }

            ReadLine();
            break;

        case E_DATA:
            if( line_.empty() || line_[0] == '#' ) { ReadLine(); return; }

            if( line_[0] == '>' ) state_ = E_READY;
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

        default: SRPRISM_ASSERT( false );
    }
}

//------------------------------------------------------------------------------
bool CFastaStream::Next(void)
{
    if( done_ ) return false;
    while( state_ != E_READY && state_ != E_EMPTY ) ProcessLine();
    if( state_ == E_EMPTY ) return false;
    seq_data_.size = seq_.size();
    state_ = E_START;
    return true;
}

END_NS( seq )
END_STD_SCOPES

