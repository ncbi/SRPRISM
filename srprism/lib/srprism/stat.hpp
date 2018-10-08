/*  $Id: stat.hpp 205525 2010-09-20 16:47:02Z morgulis $
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
 * File Description: debugging statistics
 *
 */

#ifndef __STAT_HPP__
#define __STAT_HPP__

#include <string>
#include <map>

#include "../common/def.h"
#include "../common/exception.hpp"

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CStatMap
{
    private:

        typedef std::map< std::string, size_t > TCounters;

        TCounters counters_;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode NM_COLL = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == NM_COLL ) return "name collision";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        void NewCounter( const std::string & name )
        {
            if( counters_.find( name ) == counters_.end() ) {
                counters_.insert( std::make_pair( name, (size_t)0 ) );
                return;
            }

            M_THROW( CException, NM_COLL, "counter redifinition: " << name );
        }

        size_t * GetCounter( const std::string & name )
        {
           TCounters::iterator i( counters_.find( name ) );
           if( i != counters_.end() ) return &(i->second);
           M_THROW( CException, NM_COLL, "unknown counter: " << name );
        }
};

END_NS( srprism )
END_STD_SCOPES

#endif

