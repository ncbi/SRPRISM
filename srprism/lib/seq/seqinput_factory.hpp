/*  $Id: seqinput_factory.hpp 351958 2012-02-02 15:03:39Z morgulis $
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
 * File Description: factory class for multi-column input objects
 *
 */

#ifndef __AM_SEQ_SEQINPUT_FACTORY_HPP__
#define __AM_SEQ_SEQINPUT_FACTORY_HPP__

#include "../common/def.h"

#include <string>
#include <memory>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <seq/seqinput.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>

#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CSeqInputFactory
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode TYPE = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == TYPE ) return "bad input type";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        static std::auto_ptr< CSeqInput > MakeSeqInput(
                const std::string & type,
                const std::string & name,
                int max_col,
                common::CFileBase::TCompression c );
};

END_NS( seq )
END_STD_SCOPES

#endif

