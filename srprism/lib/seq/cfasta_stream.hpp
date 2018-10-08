/*  $Id: cfasta_stream.hpp 205414 2010-09-17 17:59:42Z morgulis $
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
 * File Description: color fasta parser implementation
 *
 */

#ifndef __SEQ_CFASTA_STREAM_HPP__
#define __SEQ_CFASTA_STREAM_HPP__

#include "../common/def.h"

#include <string>

#include "../common/textfile.hpp"
#include "seqdef.hpp"
#include "stream_base.hpp"

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CColorFastaStream : public CStreamBase
{
    public:

        CColorFastaStream( 
                const std::string & name, common::CFileBase::TCompression c );

        virtual bool Next(void);

    private:

        static const TCoding INTERNAL_CODING = CODING_NCBI2NA;

        CColorFastaStream( const CColorFastaStream & );
        CColorFastaStream & operator=( const CColorFastaStream & );

        void ReadLine(void);
        void ProcessLine(void);
        void ProcessDefline(void);
        std::string Color2NCBI2NA( const std::string & cstr );

        enum EState 
        {
            E_START,
            E_DATA_START,
            E_DATA_CONT,
            E_READY
        } state_;

        TLetter lead_;
};

END_NS( seq )
END_STD_SCOPES

#endif

