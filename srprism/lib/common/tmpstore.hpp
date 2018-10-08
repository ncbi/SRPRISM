/*  $Id: tmpstore.hpp 214315 2010-12-02 21:24:25Z morgulis $
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
 * File Description: utilities for temporary file name management
 *
 */

#ifndef __AM_COMMON_TMPSTORE_HPP__
#define __AM_COMMON_TMPSTORE_HPP__

#include "../common/def.h"

/*
#ifndef NCBI_CPP_TK
#include <common/def.h>
#else
#include <../src/internal/align_toolbox/srprism/lib/common/def.h>
#endif
*/

#include <string>
#include <set>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CTmpStore
{
    public:

        CTmpStore( const std::string & tmp_dir_name = "/tmp" );
        ~CTmpStore(void);
        const std::string Register( const std::string & name );
        bool Find( const std::string & name ) const;

    private:

        CTmpStore( const CTmpStore & );
        CTmpStore & operator=( const CTmpStore & );

        typedef std::set< std::string > TData;

        const std::string CreateName( const std::string & name );

        TData data_;
        std::string tmp_name_prefix_;
        std::string tmp_name_suffix_;
};

END_NS( common )
END_STD_SCOPES

#endif

