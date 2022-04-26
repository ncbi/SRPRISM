/*  $Id: stream_factory.hpp 637057 2021-09-05 23:00:51Z morgulis $
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

#ifndef __AM_SEQ_STREAM_FACTORY_HPP__
#define __AM_SEQ_STREAM_FACTORY_HPP__

#include "../common/def.h"

#include <string>
#include <memory>

#include "../common/exception.hpp"
#include "stream_base.hpp"

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CStreamFactory
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode TYPE = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == TYPE ) return "bad input type";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        typedef int TStreamType;
        
        static const TStreamType STREAM_TYPE_FASTA   = 0;
        static const TStreamType STREAM_TYPE_FASTQ   = 1;
        static const TStreamType STREAM_TYPE_CFASTA  = 2;
        static const TStreamType STREAM_TYPE_CFASTQ  = 3;
        static const TStreamType STREAM_TYPE_ILLEGAL = 4;

        static const char * STREAM_TYPE_FASTA_NAME;
        static const char * STREAM_TYPE_FASTQ_NAME;
        static const char * STREAM_TYPE_CFASTA_NAME;
        static const char * STREAM_TYPE_CFASTQ_NAME;

        static const char * StreamTypeName( TStreamType stream_type );
        static TStreamType StreamType( const std::string & stream_type_name );

        // static std::auto_ptr< CStreamBase > MakeSeqStream(
        static std::unique_ptr< CStreamBase > MakeSeqStream(
                const std::string & stream_type_name,
                const std::string & stream_name,
                common::CFileBase::TCompression c )
        { 
            return MakeSeqStream( 
                    StreamType( stream_type_name ), stream_name, c );
        }

        // static std::auto_ptr< CStreamBase > MakeSeqStream(
        static std::unique_ptr< CStreamBase > MakeSeqStream(
                TStreamType stream_type, 
                const std::string & stream_name,
                common::CFileBase::TCompression c );

    private:

        static const char * TYPE_NAMES[];
};

END_NS( seq )
END_STD_SCOPES
#endif

