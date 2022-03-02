/*  $Id: trace.hpp 351764 2012-02-01 14:07:34Z morgulis $
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

#ifndef __AM_COMMON_TRACE_HPP__
#define __AM_COMMON_TRACE_HPP__

/*
#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif
*/

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

#include <math.h>

#include "../common/exception.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <ostream>
#include <sstream>
#include <string>

#include <math.h>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>

#endif

START_STD_SCOPES
START_NS( common )

#ifdef NDEBUG
#   define DBG_TRACE(msg)
#else
#   define DBG_TRACE(msg) M_TRACE(CTracer::DBG_LVL,msg)
#endif

#ifndef NCBI_CPP_TK

#   define M_TRACE(lvl,msg) {\
        std::ostringstream os; \
        os << msg; \
        common::CTracer::trace( __FILE__, __LINE__, lvl, os.str() ); }

#else

#   define M_TRACE(lvl,msg) {\
            switch( lvl ) {\
                case common::CTracer::DBG_LVL:     ERR_POST( Trace << msg );   break; \
                case common::CTracer::INFO_LVL:    ERR_POST( Info << msg ) ;   break; \
                case common::CTracer::WARNING_LVL: ERR_POST( Warning << msg ); break; \
                case common::CTracer::ERROR_LVL:   ERR_POST( Error << msg );   break; \
            }\
        }
#endif

class CTracer
{
    public:

        typedef int TLevel;

        static const TLevel DBG_LVL     = 0;
        static const TLevel INFO_LVL    = 1;
        static const TLevel WARNING_LVL = 2;
        static const TLevel ERROR_LVL   = 3;
        static const TLevel QUIET_LVL   = 4;

#ifndef NCBI_CPP_TK

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode NOTSET = 0;
            static const TErrorCode STREAM = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == NOTSET ) return "value not set";
                else if( code == STREAM ) return "stream error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        static void SetOutputStream( std::ostream & out );
        static void SetOutputFile( const std::string & fname );
        static void SetLevel( TLevel new_lvl );

        static std::ostream & Output( void );
        static TLevel Level( void );

        static void trace( 
                const char * file, 
                int line, 
                TLevel lvl, 
                const std::string & msg )
        { 
            if( lvl >= Curr_Lvl_ && Tr_Stream_ != 0 ) {
                (*Tr_Stream_) << Level2Str[lvl] << msg
                              << " <" << file << ":" << line << ">" 
                              << std::endl;
            }
        }

    private:

        static const TLevel NUM_LEVELS = QUIET_LVL;
        static const char * Level2Str[NUM_LEVELS];

        static std::ostream * Tr_Stream_;
        static TLevel Curr_Lvl_;
        static bool alloc_;

#endif

};

class CProgress
{
    public:

#ifndef NCBI_CPP_TK

        CProgress( 
                Uint8 total, const std::string & msg,
                double step = 1.0,
                CTracer::TLevel lvl = CTracer::ERROR_LVL ) 
            : total_( total ), current_( 0.0 ), 
              step_( step ), lvl_( lvl ), msg_( msg )
        {
            if( lvl_ >= CTracer::Level() ) {
                std::cerr << "\r" << msg_ << "0%          " << std::flush;
            }
        }

        void Current( double current )
        {
            if( lvl_ >= CTracer::Level() ) {
                double pold = 100.0*current_/total_;
                double pnew = 100.0*current/total_;

                if( fabs( pold - pnew ) >= step_ ) {
                    current_ = current;
                    std::cerr << "\r" << msg_ << std::setprecision( 4 ) 
                              << pnew << "%        " << std::flush;
                }
            }
        }

        ~CProgress()
        {
            if( lvl_ >= CTracer::Level() ) {
                std::cerr << "\r" << msg_ << "100%        " << std::endl;
            }
        }

#else

        CProgress( 
                Uint8 total, const std::string & msg,
                double step = 1.0,
                CTracer::TLevel lvl = CTracer::ERROR_LVL ) 
            : total_( total ), current_( 0.0 ), 
              step_( step ), lvl_( lvl ), msg_( msg )
        { NcbiCerr << "\r" << msg_ << "0%          " << std::flush; }

        void Current( double current )
        {
            double pold = 100.0*current_/total_;
            double pnew = 100.0*current/total_;

            if( fabs( pold - pnew ) >= step_ ) {
                current_ = current;
                NcbiCerr << "\r" << msg_ << std::setprecision( 4 ) 
                         << pnew << "%         " << std::flush;
            }
        }

        ~CProgress() 
        { NcbiCerr << "\r" << msg_ << "100%         " << std::endl; }

#endif

    private:

        CProgress( const CProgress & );
        CProgress & operator=( const CProgress & );

        Uint8 total_;
        double current_;
        double step_;
        CTracer::TLevel lvl_;
        std::string msg_;
};

class CTicker
{
    public:

#ifndef NCBI_CPP_TK

        CTicker( Uint8 step, const std::string & start_msg,
                 const std::string & end_msg = "done",
                 CTracer::TLevel lvl = CTracer::ERROR_LVL )
            : step_( step ), current_( 0 ), lvl_( lvl ), end_msg_( end_msg )
        {
            if( lvl_ >= CTracer::Level() ) {
                std::cerr << start_msg << "..." << std::flush;
            }
        }

        void Tick(void)
        {
            ++current_;
            
            if( current_ == step_ ) {
                if( lvl_ >= CTracer::Level() ) { 
                    std::cerr << '.' << std::flush; 
                }

                current_ = 0;
            }
        }

        ~CTicker(void) 
        {
            if( lvl_ >= CTracer::Level() ) {
                std::cerr << end_msg_ << std::endl;
            }
        }

#else

        CTicker( Uint8 step, const std::string & start_msg,
                 const std::string & end_msg = "done",
                 CTracer::TLevel lvl = CTracer::ERROR_LVL )
            : step_( step ), current_( 0 ), lvl_( lvl ), end_msg_( end_msg )
        { NcbiCerr << start_msg << "..." << std::flush; }

        void Tick( void ) 
        { 
            ++current_;

            if( current_ == step_ ) {
                NcbiCerr << "." << std::flush;
                current_ = 0;
            }
        }

        ~CTicker() { NcbiCerr << end_msg_ << std::endl; }

#endif

    private:

        Uint8 step_, current_;
        CTracer::TLevel lvl_;
        std::string end_msg_;
};

END_NS( common )
END_STD_SCOPES

#endif
