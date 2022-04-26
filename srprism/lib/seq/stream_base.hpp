/*  $Id: stream_base.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: base functionality for single column sequence 
 *                   input streams
 *
 */

#ifndef __SEQ_STREAM_BASE_HPP__
#define __SEQ_STREAM_BASE_HPP__

#include "../common/def.h"

#include <string>

#ifndef NCBI_CPP_TK

#include <common/textfile.hpp>
#include <seq/seqdef.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/textfile.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>

#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CStreamBase
{
    protected:

        static const TCoding CODING = CODING_IUPACNA;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode PARSE = 0;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == PARSE ) return "parse error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        typedef SSeqData< CODING, common::Uint1 > TData;

    protected:

        typedef TData::TWord TWord;
        typedef TData::TSeq  TSeq;

    public:

        typedef std::string TQual;

        CStreamBase( 
                const std::string & name, 
                common::CFileBase::TCompression c )
            : name_( name ), 
              in_( common::CReadTextFile::MakeReadTextFile( name, c ) ), 
              seq_data_( seq_, 0 ), seq_qual_( "*" ), done_( false )
        {}

        virtual ~CStreamBase() {}

        virtual bool Next( void ) = 0;
        bool Done( void ) { return done_; }

        const TSeqId & Id(void) const { return id_; }
        const TSeqTitle & Title(void) const { return title_; }
        const TData & Data(void) const { return seq_data_; }
        const TQual & Qual( void ) const { return seq_qual_; }

    protected:

        std::string name_;
        TSeq seq_;
        // std::auto_ptr< common::CReadTextFile > in_;
        std::unique_ptr< common::CReadTextFile > in_;
        TSeqId id_;
        TSeqTitle title_;
        TData seq_data_;
        TQual seq_qual_;
        std::string line_;
        bool done_;

    private:

        CStreamBase( const CStreamBase & );
        CStreamBase & operator=( const CStreamBase & );
};

END_NS( seq )
END_STD_SCOPES

#endif

