/*  $Id: textfile.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: utilities for text based file i/o
 *
 */

#ifndef __AM_COMMON_TEXTFILE_HPP__
#define __AM_COMMON_TEXTFILE_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <string>
#include <iostream>
#include <sstream>
#include <memory>

#include <common/file.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <string>
#include <iostream>
#include <sstream>
#include <memory>

#include <../src/internal/align_toolbox/srprism/lib/common/file.hpp>

#endif

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CReadTextFile : public CFileBase
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode NO_EOL = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const 
            {
                if( code == NO_EOL ) return "missing end of line";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        typedef common::Uint8 TSize;

        // static std::auto_ptr< CReadTextFile > MakeReadTextFile( 
        static std::unique_ptr< CReadTextFile > MakeReadTextFile( 
                const std::string & name, TCompression c = COMPRESSION_AUTO );

        CReadTextFile( const std::string & name )
            : CFileBase( name ), lineno_( 0 )
        {}

        virtual ~CReadTextFile() {}

        const std::string GetLine() { ++lineno_; return GetLine_Impl(); }
        virtual bool Eof() const = 0;
        TSize LineNo() const { return lineno_; }

    protected:

        virtual const std::string GetLine_Impl( void ) = 0;

        TSize lineno_;

    private:

        CReadTextFile( const CReadTextFile & );
        CReadTextFile & operator=( const CReadTextFile & );
};

class CReadTextFile_CPPStream : public CReadTextFile
{
    public:

        CReadTextFile_CPPStream( const std::string & name );
        virtual bool Eof() const { return is_.eof(); }

    protected:

        virtual const std::string GetLine_Impl( void );

    private:

        std::istream & is_;
        std::unique_ptr< std::istream > is_holder_;
        // std::auto_ptr< std::istream > is_holder_;
};


//------------------------------------------------------------------------------
class CWriteTextFile : public CFileBase
{
    public:

        // static std::auto_ptr< CWriteTextFile > MakeWriteTextFile( 
        static std::unique_ptr< CWriteTextFile > MakeWriteTextFile(
                const std::string & name, TCompression c = COMPRESSION_AUTO );

        CWriteTextFile( const std::string & name ) 
            : CFileBase( name ), lines_out_( 0 )
        {}
        
        template< typename val_t > CWriteTextFile & Out( const val_t & val )
        {
            std::ostringstream os;
            os << val;
            DoOutput( os.str() );
            return *this;
        }

        template< typename val_t > 
        CWriteTextFile & LineOut( const val_t & val )
        { Out( val ); Out( "\n" ); ++lines_out_; return *this; }

    protected:

        virtual void DoOutput( const std::string & out_str ) = 0;

        size_t lines_out_;
};

class CWriteTextFile_CPPStream : public CWriteTextFile
{
    public:

        CWriteTextFile_CPPStream( const std::string & name );

        virtual void DoOutput( const std::string & out_str )
        {
            os_ << out_str;

            if( !os_.good() ) {
                M_THROW( CFileBase::CException, WRITE,
                        "for " << name_ << " after " << lines_out_ << " lines" );
            }
        }

    private:

        std::ostream & os_;
        std::unique_ptr< std::ostream > os_holder_;
        // std::auto_ptr< std::ostream > os_holder_;
};

END_NS( common )
END_STD_SCOPES

#endif

