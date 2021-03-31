/*  $Id: util.cpp 205411 2010-09-17 17:42:11Z morgulis $
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
 * File Description: miscelanious utilities
 *
 */

#include <ncbi_pch.hpp>

#include "def.h"

#include "util.hpp"

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
const SEndianness::TTYPE SEndianness::LITTLE = 1;
const SEndianness::TTYPE SEndianness::BIG    = 0;

//------------------------------------------------------------------------------
SEndianness::TTYPE Endianness()
{
    Uint4 i = 1;
    return (((char *)&i)[0] == 0) ? SEndianness::BIG : SEndianness::LITTLE;
}

END_NS( common )
END_STD_SCOPES

