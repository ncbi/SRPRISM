/*  $Id: file.hpp 205423 2010-09-17 19:16:42Z morgulis $
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
 * File Description: utilities for file based io
 *
 */

#ifndef __AM_COMMON_FILE_HPP__
#define __AM_COMMON_FILE_HPP__

/*
#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif
*/

#ifndef NCBI_CPP_TK

#include "def.h"
#include "exception.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>
#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>

#endif

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CFileBase
{
    public:

        typedef int TCompression;

        static const TCompression COMPRESSION_AUTO = 0;
        static const TCompression COMPRESSION_NONE = 1;
        static const TCompression COMPRESSION_ZIP  = 2;
        static const TCompression COMPRESSION_BZIP = 3;

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode SYSTEM = 0;
            static const TErrorCode OPEN   = 1;
            static const TErrorCode READ   = 2;
            static const TErrorCode WRITE  = 3;
            static const TErrorCode SIZE   = 4;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == SYSTEM )     return "system error";
                else if( code == OPEN )  return "open failed";
                else if( code == READ )  return "read failed";
                else if( code == WRITE ) return "write failed";
                else if( code == SIZE )  return "length error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        typedef std::streamsize TSize;

        virtual ~CFileBase() {}

    protected:

        static TCompression Name2CompressionType( const std::string & name );

        CFileBase( const std::string & name ) : name_( name ), pos_( 0 ) {}

        const std::string name_;
        TSize pos_;

    private:

        CFileBase( const CFileBase & );
};

END_NS( common )
END_STD_SCOPES

#endif

