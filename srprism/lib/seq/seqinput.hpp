/*  $Id: seqinput.hpp 351328 2012-01-27 16:17:31Z morgulis $
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
 * File Description: base class for input parsers
 *
 */

#ifndef __AM_SEQ_SEQINPUT_HPP__
#define __AM_SEQ_SEQINPUT_HPP__

#include "../common/def.h"

#include <string>

#ifndef NCBI_CPP_TK

#include <seq/seqdef.hpp>
#include <seq/stream_base.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/stream_base.hpp>

#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CSeqInput
{
    protected:

        static const TCoding CODING = CODING_IUPACNA;

    public:

        typedef SSeqData< CODING, common::Uint1 > TData;
        typedef CStreamBase::TQual TQual;

        bool Done( void ) const { return done_; }

        virtual ~CSeqInput() {}
        virtual int NCols( void ) const = 0;
        virtual bool Next( void ) = 0;

        virtual const TData & Data( int col ) const = 0;
        virtual const TQual & Qual( int col ) const = 0;

        const TSeqId & Id( void ) const { return id_; }
        const TSeqTitle & Title( void ) const { return title_; }

    protected:

        typedef TData::TWord TWord;
        typedef TData::TSeq  TSeq;

        bool done_;
        TSeqId id_;
        TSeqTitle title_;
};

END_NS( seq )
END_STD_SCOPES

#endif

