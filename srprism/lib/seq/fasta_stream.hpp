/*  $Id: fasta_stream.hpp 205414 2010-09-17 17:59:42Z morgulis $
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
 * File Description: fasta parser implementation
 *
 */

#ifndef __SEQ_FASTA_STREAM_HPP__
#define __SEQ_FASTA_STREAM_HPP__

#include "../common/def.h"

#include <string>

#include "../common/textfile.hpp"
#include "seqdef.hpp"
#include "stream_base.hpp"

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CFastaStream : public CStreamBase
{
    public:

        CFastaStream( 
                const std::string & name, common::CFileBase::TCompression c );

        virtual bool Next(void);

    private:

        CFastaStream( const CFastaStream & );
        CFastaStream & operator=( const CFastaStream & );

        bool ReadLine(void);
        void ProcessLine(void);
        void ProcessDefline(void);

        enum EState 
        {
            E_START,
            E_DATA,
            E_READY,
            E_EMPTY
        } state_;
};

END_NS( seq )
END_STD_SCOPES

#endif

