/*  $Id: seqinput_factory.cpp 637057 2021-09-05 23:00:51Z morgulis $
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

#include <ncbi_pch.hpp>

#include "paired_stream.hpp"
#include "serial_stream.hpp"
#include "seqinput_multistream.hpp"
#include "seqinput_sam.hpp"

#ifdef USE_SRA
#   include "seqinput_sra.hpp"
#endif

#include "seqinput_factory.hpp"

START_STD_SCOPES
START_NS( seq )
USE_NS( common )

//------------------------------------------------------------------------------
// std::auto_ptr< CSeqInput > CSeqInputFactory::MakeSeqInput( 
std::unique_ptr< CSeqInput > CSeqInputFactory::MakeSeqInput( 
        const std::string & type, const std::string & name, int max_cols,
        CFileBase::TCompression c )
{
    int n_names( 0 );

    if( !name.empty() ) {
        std::string::size_type epos( 0 );

        while( true ) {
            ++n_names;

            if( (epos = name.find_first_of( ",", epos )) == 
                    std::string::npos ) {
                break;
            }

            ++epos;
        }
    }
    else n_names = 1;

    if( n_names > 1 && (type == "sam" || type == "sra") )
    {
        M_THROW(
            CException, FORMAT,
            "multiple inputs are not supported with input format " << type );
    }

    if( type == "sam" ) {
        // return std::auto_ptr< CSeqInput >(
        return std::unique_ptr< CSeqInput >(
                new CSeqInput_SAM( name, max_cols == 2, c ) );
    }
#ifdef USE_SRA
    else if( type == "sra" )
    {
        return std::auto_ptr< CSeqInput >(
                new CSeqInput_SRA( name, max_cols ) );
    }
#endif
    else if( n_names == 1 && max_cols == 2 ) {
        return MakePairedStream( type, name, c );
    }
    else if( n_names > 1 && max_cols == 1 )
    {
        return MakeSerialStream( type, name, c );
    }
    // else return std::auto_ptr< CSeqInput >( 
    else return std::unique_ptr< CSeqInput >( 
            new CSeqInputMultiStream( name, type, max_cols, c ) );
}

END_NS( seq )
END_STD_SCOPES

