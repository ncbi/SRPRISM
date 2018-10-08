/*  $Id: text_formatter.hpp 205411 2010-09-17 17:42:11Z morgulis $
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
 * File Description: simple text formatter to format usage text
 *
 */

#ifndef __AM_COMMON_TEXT_FORMATTER_HPP_
#define __AM_COMMON_TEXT_FORMATTER_HPP_

#include "def.h"

#include <string>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class CTextFormatter
{
    public:

        CTextFormatter( Uint1 start, Uint1 len )
            : prefix_( start, ' ' ), len_( len )
        {}

        const std::string operator()( const std::string & text ) const;

    private:

        const std::string prefix_;
        Uint1 len_;
};

END_NS( common )
END_STD_SCOPES

#endif

