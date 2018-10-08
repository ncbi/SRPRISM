/*  $Id: text_formatter.cpp 205411 2010-09-17 17:42:11Z morgulis $
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

#include <ncbi_pch.hpp>

#include "def.h"
#include "text_formatter.hpp"

#include <iostream>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
static inline void _add_word( 
        std::string & result, std::string & word, Uint1 & c )
{ result.append( word ); c += word.size(); word.clear(); }
                               
//------------------------------------------------------------------------------
const std::string CTextFormatter::operator()( const std::string & text ) const
{
    typedef std::string::const_iterator TIter;
    Uint1 c = 0;
    std::string result;
    std::string word;
    std::string curr_prefix = prefix_;
    bool line_start = true;

    for( TIter i_ltr = text.begin(); i_ltr != text.end(); ++i_ltr ) {
        if( *i_ltr == '\t' ) { curr_prefix.append( 4, ' ' ); continue; }
        if( line_start ) { result.append( curr_prefix ); line_start = false; }

        switch( *i_ltr ) {
            case '\n':
                
                if( c + word.size() > len_ ) {
                    result.push_back( '\n' ); c = 0;
                    result.append( curr_prefix );
                }

                _add_word( result, word, c );
                result.append( "\n\n" );
                curr_prefix = prefix_;
                line_start = true;
                c = 0;
                break;

            case ' ': 

                if( c + word.size() > len_ ) { 
                    result.push_back( '\n' ); c = 0; 
                    result.append( curr_prefix ); 
                }

                _add_word( result, word, c );
                result.push_back( ' ' );
                ++c;
                break;

            default: word.push_back( *i_ltr ); break;
        }
    }

    if( c + word.size() > len_ ) { 
        result.push_back( '\n' );
        result.append( curr_prefix );
    }

    result.append( word );
    return result;
}

END_NS( common )
END_STD_SCOPES

