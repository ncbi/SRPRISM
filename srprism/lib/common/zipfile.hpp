/*  $Id: zipfile.hpp 205312 2010-09-16 19:22:21Z morgulis $
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
 * File Description: support for zlib based text file compression
 *
 */

#ifndef __AM_COMMON_ZIPFILE_HPP__
#define __AM_COMMON_ZIPFILE_HPP__

#include <common/def.h>

#ifndef WIN32
#include <zlib.h>
#endif

#include <common/exception.hpp>
#include <common/textfile.hpp>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CReadTextFile_Zip : public CReadTextFile
{
    public:

        struct CException : public CReadTextFile::CException
        {
            typedef CReadTextFile::CException TBase;

            static const TErrorCode ZIP_ERROR = 1;
            static const TErrorCode IO_ERROR  = 2;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                switch( code ) {
                    case ZIP_ERROR: return "zlib error";
                    case IO_ERROR: return "system i/o error";
                    default: return TBase::ErrorMessage( code );
                }
            }

            M_EXCEPT_CTOR_2( CException, TBase );
        };

        CReadTextFile_Zip( const std::string & name );

#ifndef WIN32
        virtual ~CReadTextFile_Zip() { gzclose( gzf_ ); }
        virtual bool Eof() const { return (gzeof( gzf_ ) == 1); }
#else
        virtual ~CReadTextFile_Zip() {}
        virtual bool Eof() const { return true; }
#endif

    protected:

        virtual const std::string GetLine_Impl( void );

    private:

#ifndef WIN32
        gzFile gzf_;
#endif
};

END_NS( common )
END_STD_SCOPES

#endif

