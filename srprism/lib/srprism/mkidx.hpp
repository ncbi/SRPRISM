/*  $Id: mkidx.hpp 384196 2012-12-21 17:04:48Z morgulis $
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
 * File Description: index creation code
 *
 */

#ifndef __SRPRIZM_MKIDX_HPP__
#define __SRPRIZM_MKIDX_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <string>
#include <vector>

#include <common/exception.hpp>
#include <common/binfile.hpp>
#include <srprism/index_base.hpp>

#else

#include <string>

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/binfile.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/index_base.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CMkIdx : public CIndexBase
{
    public:

        struct SOptions
        {
            SOptions()
                : infmt( "fasta" ), outfmt( "standard" ),
                  input_compression( common::CFileBase::COMPRESSION_AUTO ),
                  max_mem( 2048 ), ss_seg_len( 8192 ), al_extend( 2000 )
            {
            }

            std::string input;
            std::string input_list;
            std::string alt_loc_spec;
            std::string output;
            std::string infmt;
            std::string outfmt;
            common::CFileBase::TCompression input_compression;
            size_t max_mem;
            common::Uint4 ss_seg_len;
            size_t al_extend;
        };

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode VALIDATE = 0;
            static const TErrorCode MEMORY   = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == VALIDATE ) return "validation error";
                else if( code == MEMORY ) return "out of memory";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CMkIdx( const SOptions & options );
        
        void Run( void );

    private:

        static void SaveIdxHeader( common::CWriteBinFile & idx_file );
        void Validate( void );
        void MkSeqStore( void );

        std::vector< std::string > input_;
        std::string alt_loc_spec_name_;
        std::string output_;
        std::string infmt_;
        std::string outfmt_;
        common::CFileBase::TCompression input_c_;
        size_t max_mem_;
        size_t ss_seg_len_;
        size_t al_extend_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

