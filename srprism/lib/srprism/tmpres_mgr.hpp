/*  $Id: tmpres_mgr.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: temporary result storage manipulation
 *
 */

#ifndef __SRPRISM_TMPRES_MGR_HPP__
#define __SRPRISM_TMPRES_MGR_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <cstdlib>
#include <memory>

#include <common/exception.hpp>
#include <common/tmpstore.hpp>
#include <common/binfile.hpp>
#include <srprism/result.hpp>

#else

#include <cstdlib>
#include <memory>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/tmpstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/binfile.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CTmpResMgr
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode READ = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == READ ) return "read error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CTmpResMgr( 
                void * mainbuf, size_t mainbuf_size, 
                const std::string & tmp_name, common::CTmpStore & tmp_store );

        CResult Save( size_t res_len )
        {
            if( mainbuf_.Full( res_len ) ) {
                CheckWriteTmpFile();
                mainbuf_.Dump( *os_ );
            }

            return mainbuf_.Add( res_len );
        }

        void LoadInit(void);
        void LoadFinal(void);
        CResult Load( void );

    private:

        class CTmpResBuf
        {
            public:

                CTmpResBuf( void * buf, size_t buf_size );

                bool Full( size_t sz ) const 
                { return (curr_size_ + sz >= bufsize_); }

                bool Last(void) const;

                CResult Add( size_t res_len ) 
                { 
                    CResult r( buf_ + curr_size_ );
                    curr_size_ += res_len; 
                    return r;
                }

                void Dump( common::CWriteBinFile & os );

                void ReadInit(void) { read_pos_ = 0; }
                void WriteInit(void) { curr_size_ = 0; }

                bool Load( common::CReadBinFile & is );
                CResult ReadNext( void );

            private:

                CTmpResBuf( const CTmpResBuf & );
                CTmpResBuf & operator=( const CTmpResBuf & );

                char * buf_;
                size_t bufsize_, curr_size_, read_pos_;
        };

        CTmpResMgr( const CTmpResMgr & );
        CTmpResMgr & operator=( const CTmpResMgr & );

        void CheckWriteTmpFile(void);
        bool CheckReadTmpFile(void);

        common::CTmpStore & tmp_store_;
        CTmpResBuf mainbuf_;
        const std::string tmp_name_;
        std::unique_ptr< common::CWriteBinFile > os_;
        std::unique_ptr< common::CReadBinFile > is_;

        /*
        std::auto_ptr< common::CWriteBinFile > os_;
        std::auto_ptr< common::CReadBinFile > is_;
        */
};

END_NS( srprism )
END_STD_SCOPES

#endif

