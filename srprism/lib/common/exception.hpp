/*  $Id: exception.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: exception management
 *
 */

#ifndef __AM_COMMON_EXCEPTION_HPP__
#define __AM_COMMON_EXCEPTION_HPP__

#include "../common/def.h"

/*
#ifndef NCBI_CPP_TK

#include <common/def.h>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#endif
*/

#include <string>
#include <sstream>
#include <stdexcept>

START_STD_SCOPES
START_NS( common )

#ifndef NCBI_CPP_TK

#ifndef M_THROW
#   define M_THROW(eclass,code,msg) { \
        std::ostringstream os; \
        os << msg << std::flush; \
        throw eclass( __FILE__, __LINE__, eclass::code, os.str() ); \
    }
#endif

#ifndef M_EXCEPT_CTOR
#   define M_EXCEPT_CTOR(eclass) eclass( \
        const std::string & file, \
        unsigned long line, \
        TErrorCode err_code, \
        const std::string & msg ) \
        : common::CException( file, line, err_code, msg ) { \
            msg_ += ErrorMessage( code_) + " (" + msg + ")"; }
#endif

#ifndef M_EXCEPT_CTOR_2
#   define M_EXCEPT_CTOR_2(eclass,bclass) eclass( \
        const std::string & file, \
        unsigned long line, \
        TErrorCode err_code, \
        const std::string & msg ) \
        : bclass( file, line, err_code, msg ) { \
            msg_ += ErrorMessage( code_) + " (" + msg + ")"; }
#endif

//------------------------------------------------------------------------------

class CException : public std::runtime_error
{
    public:

        typedef Uint4 TErrorCode;

        static const TErrorCode NOCODE = 0;

        CException( const std::string & file,
                    unsigned long line,
                    TErrorCode err_code,
                    const std::string & msg )
            : std::runtime_error( "" ), code_( err_code )
        {
            std::ostringstream os;
            os << file << ":" << line << " [" << code_ << "] ";
            msg_ = os.str();
        }

        virtual ~CException() throw() {}

        virtual const std::string ErrorMessage( TErrorCode ) const
        { return "unspecified exception"; }

        const char * what() const throw() { return msg_.c_str(); }

        TErrorCode ErrorCode() const { return code_; }

    protected:

        TErrorCode code_;
        std::string msg_;
};

#else

#include <corelib/ncbiexpt.hpp>

#ifndef M_THROW
#   define M_THROW(eclass,code,msg) { \
        std::ostringstream os; \
        os << msg << std::flush; \
        NCBI_THROW( eclass, code, os.str() ); \
    }
#endif

#ifndef M_EXCEPT_CTOR
#   define M_EXCEPT_CTOR(eclass) eclass( \
        ncbi::CDiagCompileInfo inf, \
        const ncbi::CException * pe, \
        TErrorCode err_code, \
        const std::string & msg ) \
        : common::CException( inf, pe, (ncbi::CException::EErrCode)err_code, msg ) \
        { base_msg_ = ErrorMessage( err_code ); \
          base_msg_ += ": " + msg; }
#endif

#ifndef M_EXCEPT_CTOR_2
#   define M_EXCEPT_CTOR_2(eclass,bclass) eclass( \
        ncbi::CDiagCompileInfo inf, \
        const ncbi::CException * pe, \
        TErrorCode err_code, \
        const std::string & msg ) \
        : bclass( inf, pe, (ncbi::CException::EErrCode)err_code, msg ) \
        { base_msg_ = ErrorMessage( err_code ); \
          base_msg_ += ": " + msg; }
#endif

class CException : public ncbi::CException
{
    public:

        typedef Uint4 TErrorCode;
        typedef ncbi::CException TBase;

        static const TErrorCode NOCODE = 0;

        CException( const std::string & file,
                    unsigned long line,
                    TErrorCode err_code,
                    const std::string & msg )
            : TBase( ncbi::CDiagCompileInfo( file.c_str(), line ), 0,
                     (TBase::EErrCode)err_code, msg ), code_( err_code )
        {}

        CException( ncbi::CDiagCompileInfo inf, const TBase * pe, 
                    TBase::EErrCode err_code, const string & msg )
            : TBase( inf, pe, err_code, msg ), code_( err_code )
        {}

        virtual ~CException() throw() {}

        virtual const std::string ErrorMessage( TErrorCode ) const
        { return "unspecified exception"; }

        TErrorCode ErrorCode() const { return code_; }

        virtual const char * GetErrCodeString( void ) const
        { return base_msg_.c_str(); }

        const char * what() const throw() { return base_msg_.c_str(); }

    protected:

        TErrorCode code_;
        std::string base_msg_;
};

#endif

END_NS( common )
END_STD_SCOPES

#endif

