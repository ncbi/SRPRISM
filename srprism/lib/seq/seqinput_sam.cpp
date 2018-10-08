/*  $Id: seqinput_sam.cpp 351958 2012-02-02 15:03:39Z morgulis $
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
 * File Description: input parser for sam-like inputs
 *
 */

#include <ncbi_pch.hpp>

#include "seqinput_sam.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
CSeqInput_SAM::CSeqInput_SAM( 
        const std::string & name, bool paired, 
        common::CFileBase::TCompression c )
    : paired_( paired ), name_( name ),
      in_( CReadTextFile::MakeReadTextFile( name, c ) ),
      d_0_( s_[0], 0 ), d_1_( s_[1], 0 )
{
    done_ = false;
    SkipHeader();
}

//------------------------------------------------------------------------------
void CSeqInput_SAM::ReadLine( void )
{
    if( in_->Eof() ) { line_ = ""; return; }
    line_ = in_->GetLine();

    if( !in_->Eof() && line_.empty() ) {
        M_THROW( CException, PARSE, 
                 "at line " << in_->LineNo() << 
                 ": empty lines are not supported" );
    }
}

//------------------------------------------------------------------------------
void CSeqInput_SAM::SkipHeader( void )
{
    if( !in_->Eof() ) {
        do {
            ReadLine();
            if( !line_.empty() && line_[0] != '@' ) break;
        }while( !line_.empty() );
    }
}

//------------------------------------------------------------------------------
void CSeqInput_SAM::MakeWords( void )
{
    words_.clear();
    std::string::size_type spos, epos;
    spos = line_.find_first_not_of( "\t ", 0 );

    if( spos != 0 ) {
        M_THROW( CException, PARSE,
                 "at line " << in_->LineNo() << 
                 ": line can not start with space" );
    }

    while( spos != std::string::npos ) {
        epos = line_.find_first_of( "\t ", spos );
        words_.push_back( line_.substr( spos, epos - spos ) );
        spos = line_.find_first_not_of( "\t ", epos );
    }

    if( words_.size() < MIN_FIELDS ) {
        M_THROW( CException, PARSE,
                 "at line " << in_->LineNo() << 
                 ": entry contains too few fields" );
    }
}

//------------------------------------------------------------------------------
CSeqInput_SAM::SEntry CSeqInput_SAM::ParseEntry( void )
{
    MakeWords();
    SEntry res;
    res.id = words_[ID_FIELD];
    res.seq = words_[SEQ_FIELD];
    res.qual = words_[QUAL_FIELD];
    res.flag = atoi( words_[FLAG_FIELD].c_str() );

    if( res.flag == FLAG_SINGLE || res.flag == FLAG_PAIRED_1 ) res.idx = 0;
    else if( res.flag == FLAG_PAIRED_2 ) res.idx = 1;
    else M_THROW( CException, PARSE,
                  "at line " << in_->LineNo() <<
                  "only the following flag values are recognized: [" <<
                  FLAG_SINGLE << ',' << 
                  FLAG_PAIRED_1 << ',' << FLAG_PAIRED_2 << ']' );

    return res;
}

//------------------------------------------------------------------------------
std::string CSeqInput_SAM::CheckIDs( 
        const std::string & id_1, const std::string & id_2 )
{
    size_t sz( 0 );

    {
        size_t min_sz( std::min( id_1.size(), id_2.size() ) );
        for( ; sz < min_sz; ++sz ) if( id_1[sz] != id_2[sz] ) break;
    }

    bool match( false );
    
    if( sz > 0 && (
                id_1[sz-1] == '/' ||
                id_1[sz-1] == '.' ||
                id_1[sz-1] == '_' ) ) {
        match = true;

        for( size_t i( sz ); i < id_1.size(); ++i ) {
            if( !isdigit( id_1[sz] ) ) { match = false; break; }
        }

        for( size_t i( sz ); i < id_2.size(); ++i ) {
            if( !isdigit( id_2[sz] ) ) { match = false; break; }
        }

        if( match ) {
            size_t id_num( (size_t)atoi( id_1.substr( sz ).c_str() ) );
            if( id_num != 1 ) match = false;
            id_num = (size_t)atoi( id_1.substr( sz ).c_str() );
            if( id_num != 2 ) match = false;
        }

        if( match ) --sz;
    }

    if( !match ) {
        if( sz < id_1.size() || sz < id_2.size() ) {
            M_THROW( CException, PARSE,
                     "at line " << in_->LineNo() <<
                     ": mate ids do not match" );
        }

        return id_1;
    }

    if( sz == 0 ) {
        M_THROW( CException, PARSE,
                 "at line " << in_->LineNo() <<
                 ": mate ids do not match" );
    }

    return id_1.substr( 0, sz );
}

//------------------------------------------------------------------------------
bool CSeqInput_SAM::NextPaired( void )
{
    if( !line_.empty() ) {
        SEntry e0( ParseEntry() );
        ReadLine();

        if( line_.empty() ) {
            M_THROW( CException, PARSE, "odd number of lines in the input" );
        }

        SEntry e1( ParseEntry() );

        if( (e0.flag != FLAG_PAIRED_1 && e0.flag != FLAG_PAIRED_2) ||
                (e1.flag != FLAG_PAIRED_1 && e1.flag != FLAG_PAIRED_2) ) {
            M_THROW( CException, PARSE, 
                     "at line " << in_->LineNo() << ": " <<
                     "wrong flag value for paired-end input; " <<
                     "allowed values are " << 
                     FLAG_PAIRED_1 << "," << FLAG_PAIRED_2 );
        }

        title_.clear();
        id_ = CheckIDs( e0.id, e1.id );
        d_0_.size = (e0.idx == 0) ? e0.seq.size() : e1.seq.size();
        d_1_.size = (e0.idx == 0) ? e1.seq.size() : e0.seq.size();
        s_[e0.idx].resize( e0.seq.size() );
        s_[e1.idx].resize( e1.seq.size() );
        std::copy( e0.seq.begin(), e0.seq.end(), s_[e0.idx].begin() );
        std::copy( e1.seq.begin(), e1.seq.end(), s_[e1.idx].begin() );
        q_[e0.idx] = e0.qual;
        q_[e1.idx] = e1.qual;
        ReadLine();
        return true;
    }
    else{ done_ = true; return false; }
}

//------------------------------------------------------------------------------
bool CSeqInput_SAM::NextUnpaired( void )
{
    if( !line_.empty() ) {
        SEntry e( ParseEntry() );

        if( e.flag != FLAG_SINGLE ) {
            M_THROW( CException, PARSE, 
                     "at line " << in_->LineNo() << ": " <<
                     "wrong flag value for non-paired input; " <<
                     "allowed value is " << FLAG_SINGLE );
        }

        s_[0].resize( e.seq.size() );
        std::copy( e.seq.begin(), e.seq.end(), s_[0].begin() );
        d_0_.size = e.seq.size();
        q_[0] = e.qual;
        id_ = e.id;
        title_.clear();
        ReadLine();
        return true;
    }
    else{ done_ = true; return false; }
}

//------------------------------------------------------------------------------
bool CSeqInput_SAM::Next( void )
{
    return paired_ ? NextPaired() : NextUnpaired();
}

END_NS( seq )
END_STD_SCOPES

