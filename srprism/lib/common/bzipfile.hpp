/*  $Id: bzipfile.hpp 205312 2010-09-16 19:22:21Z morgulis $
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
 * File Description: support for bzip2 based text file compression
 *
 */

#ifndef __AM_COMMON_BZIPFILE_HPP__
#define __AM_COMMON_BZIPFILE_HPP__

#include <common/def.h>

#include <cstdio>

#ifndef WIN32
#include <bzlib.h>
#endif

#include <common/exception.hpp>
#include <common/textfile.hpp>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CReadTextFile_BZip : public CReadTextFile
{
    public:

        struct CException : public CReadTextFile::CException
        {
            typedef CReadTextFile::CException TBase;

            static const TErrorCode BZIP_ERROR = 1;
            static const TErrorCode IO_ERROR   = 2;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                switch( code ) {
                    case BZIP_ERROR: return "bzlib error";
                    case IO_ERROR: return "system i/o error";
                    default: return TBase::ErrorMessage( code );
                }
            }

            M_EXCEPT_CTOR_2( CException, TBase );
        };

        CReadTextFile_BZip( const std::string & name );

        virtual ~CReadTextFile_BZip() 
        { 
#ifndef WIN32
            int err( 0 );
            BZ2_bzReadClose( &err, bzf_ );
            fclose( f_ );
#endif
        }

        virtual bool Eof() const { return eof_; }

    protected:

        virtual const std::string GetLine_Impl( void );

    private:

        bool eof_;
        FILE * f_;
#ifndef WIN32
        BZFILE * bzf_;
#endif
};

END_NS( common )
END_STD_SCOPES

#endif

