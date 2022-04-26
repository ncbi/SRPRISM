/*  $Id: stream_factory.cpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: factory class for single column input stream objects
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include "stream_factory.hpp"
#include "fasta_stream.hpp"
#include "fastq_stream.hpp"
#include "cfasta_stream.hpp"
#include "cfastq_stream.hpp"

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
const char * CStreamFactory::STREAM_TYPE_FASTA_NAME  = "fasta";
const char * CStreamFactory::STREAM_TYPE_FASTQ_NAME  = "fastq";
const char * CStreamFactory::STREAM_TYPE_CFASTA_NAME = "cfasta";
const char * CStreamFactory::STREAM_TYPE_CFASTQ_NAME = "cfastq";

const char * CStreamFactory::TYPE_NAMES[] = {
    STREAM_TYPE_FASTA_NAME,
    STREAM_TYPE_FASTQ_NAME,
    STREAM_TYPE_CFASTA_NAME,
    STREAM_TYPE_CFASTQ_NAME,
    0
};

//------------------------------------------------------------------------------
const char * CStreamFactory::StreamTypeName( TStreamType type )
{ return TYPE_NAMES[type]; }

//------------------------------------------------------------------------------
CStreamFactory::TStreamType CStreamFactory::StreamType(
        const std::string & type_name )
{
    for( int i( 0 ); i < STREAM_TYPE_ILLEGAL; ++i ) {
        if( type_name == TYPE_NAMES[i] ) return i;
    }

    return STREAM_TYPE_ILLEGAL;
}

//------------------------------------------------------------------------------
// std::auto_ptr< CStreamBase > CStreamFactory::MakeSeqStream(
std::unique_ptr< CStreamBase > CStreamFactory::MakeSeqStream(
        TStreamType type, const std::string & name, 
        common::CFileBase::TCompression c )
{
    switch( type ) {
        case STREAM_TYPE_FASTA:
            // return std::auto_ptr< CStreamBase >(
            return std::unique_ptr< CStreamBase >(
                    new CFastaStream( name, c ) );
        case STREAM_TYPE_FASTQ:
            // return std::auto_ptr< CStreamBase >(
            return std::unique_ptr< CStreamBase >(
                    new CFastqStream( name, c ) );
        case STREAM_TYPE_CFASTA:
            // return std::auto_ptr< CStreamBase >(
            return std::unique_ptr< CStreamBase >(
                    new CColorFastaStream( name, c ) );
        case STREAM_TYPE_CFASTQ:
            // return std::auto_ptr< CStreamBase >(
            return std::unique_ptr< CStreamBase >(
                    new CColorFastqStream( name, c ) );
    }

    M_THROW( CException, TYPE, "" );
}

END_NS( seq )
END_STD_SCOPES

