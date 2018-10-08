/*  $Id: search_mode.hpp 384196 2012-12-21 17:04:48Z morgulis $
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
 * File Description: search mode related parameters
 *
 */

#ifndef __SRPRISM_SEARCH_MODE_HPP__
#define __SRPRISM_SEARCH_MODE_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <srprism/scoring.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/srprism/scoring.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
struct SSearchMode
{
    static const int DEFAULT   = 0;
    static const int PARTIAL   = 1;
    static const int SUM_ERR   = 2;
    static const int BOUND_ERR = 3;
};

//------------------------------------------------------------------------------
template< int smode >
struct SSearchModeTraits
{
    typedef CDefaultScoringSystem TScoringSys;
    static const bool PARTIAL_ALIGNMENT = false;
};

template<> struct SSearchModeTraits< SSearchMode::PARTIAL >
{
    typedef CLenMinusErrScoringSystem TScoringSys;
    static const bool PARTIAL_ALIGNMENT = true;
};

template<> struct SSearchModeTraits< SSearchMode::SUM_ERR >
{
    typedef CSumErrScoringSystem TScoringSys;
    static const bool PARTIAL_ALIGNMENT = false;
};

template<> struct SSearchModeTraits< SSearchMode::BOUND_ERR >
{
    typedef CWeakScoringSystem TScoringSys;
    static const bool PARTIAL_ALIGNMENT = false;
};

END_NS( srprism )
END_STD_SCOPES

#endif

