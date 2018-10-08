/*  $Id: trace.cpp 205515 2010-09-20 15:51:04Z morgulis $
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
 * File Description: utilities for log management
 *
 */

#include <ncbi_pch.hpp>

#include "def.h"

#include <fstream>

#include "trace.hpp"

START_STD_SCOPES
START_NS( common )

#ifndef NCBI_CPP_TK

//------------------------------------------------------------------------------
const char * CTracer::Level2Str[CTracer::NUM_LEVELS] = {
    "DEBUG:   ", 
    "INFO:    ", 
    "WARNING: ", 
    "ERROR:   " };
std::ostream * CTracer::Tr_Stream_ = 0;
CTracer::TLevel CTracer::Curr_Lvl_ = CTracer::DBG_LVL;
bool CTracer::alloc_ = false;

//------------------------------------------------------------------------------
void CTracer::SetOutputStream( std::ostream & out ) 
{ 
    if( alloc_ ) delete Tr_Stream_;
    Tr_Stream_ = &out; alloc_ = false;
}

//------------------------------------------------------------------------------
void CTracer::SetOutputFile( const std::string & fname )
{
    if( alloc_ ) delete Tr_Stream_;
    Tr_Stream_ = new std::ofstream( fname.c_str() );

    if( !Tr_Stream_ ) {
        M_THROW( CException, STREAM, "failed to open " << fname );
    }

    alloc_ = true;
}

//------------------------------------------------------------------------------
void CTracer::SetLevel( TLevel new_lvl ) { Curr_Lvl_ = new_lvl; }

//------------------------------------------------------------------------------
std::ostream & CTracer::Output( void )
{
    if( Tr_Stream_ != 0 ) return *Tr_Stream_;
    else M_THROW( CException, NOTSET, "no trace stream" );
}

//------------------------------------------------------------------------------
CTracer::TLevel CTracer::Level( void ) { return Curr_Lvl_; }

#endif

END_NS( common )
END_STD_SCOPES

