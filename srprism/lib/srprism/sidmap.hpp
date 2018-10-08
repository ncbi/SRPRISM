/*  $Id: sidmap.hpp 351764 2012-02-01 14:07:34Z morgulis $
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
 * File Description: subject sequence id map interface
 *
 */

#ifndef __SRPRISM_SIDMAP_HPP__
#define __SRPRISM_SIDMAP_HPP__

#include "../common/def.h"

#include <string>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/seqstore.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CSIdMap
{
    typedef TDBOrdId TSeqId;

    static const char * FILE_SFX;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode FORMAT = 0;
            static const TErrorCode MEMORY = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == FORMAT ) return "format error";
                else if( code == MEMORY ) return "out of memory";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CSIdMap( const std::string & name, CMemoryManager & mem_mgr );
        ~CSIdMap(void);

        const std::string operator[]( TSeqId id ) const
        { 
            return std::string( 
                    data_ + offs_[id], offs_[id + 1] - offs_[id] ); 
        }

    private:

        CSIdMap( const CSIdMap & );
        CSIdMap & operator=( const CSIdMap & );

        void CleanUp(void);

        CMemoryManager & mem_mgr_;

        size_t * offs_;
        char * data_;
        size_t n_ids_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

