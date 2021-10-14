/*  $Id: seqinput_multistream.hpp 639115 2021-10-13 15:24:22Z morgulis $
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
 * File Description: input parser that reads multiple input streams
 *                   in parallel as columns
 *
 */

#ifndef __AM_SEQ_SEQINPUT_MULTISTREAM_HPP__
#define __AM_SEQ_SEQINPUT_MULTISTREAM_HPP__

#include "../common/def.h"

#include <string>
#include <vector>

#include "stream_base.hpp"
#include "seqinput.hpp"

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CSeqInputMultiStream : public CSeqInput
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode EOS = 0;
            static const TErrorCode COL = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == EOS ) return "end of stream error";
                else if( code == COL ) return "inconsistent columns";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CSeqInputMultiStream( 
                const std::string & name, 
                const std::string & type, 
                int max_cols,
                common::CFileBase::TCompression c );

        virtual ~CSeqInputMultiStream( void );
        virtual int NCols( void ) const { return (int)streams_.size(); }

        virtual const TData & Data( int col ) const 
        { return streams_[col]->Data(); }

        virtual const TQual & Qual( int col ) const
        { return streams_[col]->Qual(); }

        virtual bool Next( void );

    private:

        CSeqInputMultiStream( const CSeqInputMultiStream & );
        CSeqInputMultiStream & operator=( const CSeqInputMultiStream & );

        typedef std::vector< CStreamBase * > TStreams;

        TStreams streams_;
        common::Uint8 seq_n_;
};

END_NS( seq )
END_STD_SCOPES

#endif

