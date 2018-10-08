/*  $Id: seqinput_sra.hpp 540690 2017-07-10 15:23:27Z morgulis $
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
 */

#ifndef AM_SEQ_SEQINPUT_SRA_HPP
#define AM_SEQ_SEQINPUT_SRA_HPP

#include "../common/def.h"

#ifdef USE_SRA

#include <ngs/ReadCollection.hpp>

#ifndef NCBI_CPP_TK
#   include <seq/seqinput.hpp>
#else
#   include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#endif

START_STD_SCOPES
START_NS( seq )

using namespace ngs;

//==============================================================================
class CSeqInput_SRA : public CSeqInput
{
    static std::string const QUAL;

public:

    struct CException : public common::CException
    {
        typedef common::CException TBase;

        virtual const std::string ErrorMessage( TErrorCode code ) const
        {
            return TBase::ErrorMessage( code );
        }

        M_EXCEPT_CTOR( CException )
    };

    CSeqInput_SRA( std::string const & acc, int n_cols );
    virtual ~CSeqInput_SRA() {}

    virtual int NCols( void ) const { return n_cols_; }
    virtual bool Next( void );
    virtual const TQual & Qual( int col ) const { return QUAL; }

    virtual const TData & Data( int col ) const
    {
        return col == 0 ? d0_ : d1_;
    }

private:

    std::string acc_;
    int n_cols_,
        n_frags_,
        c_col_;
    TData d0_, d1_;
    TSeq s0_, s1_;

    ReadCollection run_;
    ReadIterator it_;

    size_t start_;
};

END_NS( seq )
END_STD_SCOPES

#endif

#endif

