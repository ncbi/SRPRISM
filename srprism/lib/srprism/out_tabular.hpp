/*  $Id: out_tabular.hpp 431273 2014-04-02 17:10:44Z morgulis $
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
 * File Description: output in tabular format
 *
 */

#ifndef __SRPRISM_OUT_TABULAR_HPP__
#define __SRPRISM_OUT_TABULAR_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

// #include <common/def.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <common/exception.hpp>

#include <srprism/out_base.hpp>
#include <srprism/result.hpp>

#else

// #include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>

#include <../src/internal/align_toolbox/srprism/lib/srprism/out_base.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class COutTabular : public COutBase
{
    static const std::string SEP;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode OPEN = 0;

            const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == OPEN ) return "open error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException );
        };

        COutTabular( 
                const std::string & name,
                const std::string & input,
                const std::string & input_fmt,
                common::CFileBase::TCompression input_c,
                bool skip_unmapped,
                bool force_paired,
                bool force_unpaired,
                bool no_qids,
                CSeqStore * seq_store,
                CSIdMap * sid_map )
            : COutBase( 
                    name, input, input_fmt, input_c,  
                    skip_unmapped, force_paired, force_unpaired, no_qids,
                    seq_store, sid_map )
        {
        }

    private:

        virtual void ResultOut( 
                const CResult & result, bool map_unmapped, 
                TQueryOrdId q_adj, int const * , bool primary = true );

        COutTabular( const COutTabular & );
        COutTabular & operator=( const COutTabular & );
};

END_NS( srprism )
END_STD_SCOPES

#endif

