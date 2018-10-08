/*  $Id: binfile.hpp 205423 2010-09-17 19:16:42Z morgulis $
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
 * File Description: utilities for binary file i/o.
 *
 */

#ifndef __AM_COMMON_BINFILE_HPP__
#define __AM_COMMON_BINFILE_HPP__

#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <string>
#include <fstream>

#include "../common/exception.hpp"
#include "../common/file.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <string>
#include <fstream>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/file.hpp>

#endif

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CReadBinFile : public CFileBase
{
    public:

        CReadBinFile( const std::string & name );

        TSize Read( char * buf, TSize n, bool strict = false );
        bool Eof() const { return is_.eof(); }

    private:

        CReadBinFile( const CReadBinFile & );
        CReadBinFile & operator=( const CReadBinFile & );

        std::ifstream is_;
};

//------------------------------------------------------------------------------
class CWriteBinFile : public CFileBase
{
    public:

        CWriteBinFile( const std::string & name );

        void Write( const char * buf, TSize n );
        size_t BytesWritten( void ) const { return pos_; }

    private:

        CWriteBinFile( const CWriteBinFile & );
        CWriteBinFile & operator=( const CWriteBinFile & );

        std::ofstream os_;
};

END_NS( common )
END_STD_SCOPES

#endif

