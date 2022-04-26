/*  $Id: seqinput_sam.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: input parser for sam-like inputs
 *
 */

#ifndef __AM_SEQ_SEQINPUT_SAM_HPP__
#define __AM_SEQ_SEQINPUT_SAM_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK
#   include <common/textfile.hpp>
#   include <seq/seqinput.hpp>
#else
#   include <../src/internal/align_toolbox/srprism/lib/common/textfile.hpp>
#   include <../src/internal/align_toolbox/srprism/lib/seq/seqinput.hpp>
#endif

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
class CSeqInput_SAM : public CSeqInput
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode PARSE = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == PARSE ) return "parse error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

    private:

        static const int MAX_COLS = 2;

#ifndef NDEBUG
        void CheckCols( int col ) const
        {
            assert( col < MAX_COLS );
            assert( paired_ || col == 0 );
        }
#   define CHECK_COLS(c) CheckCols( c )
#else
#   define CHECK_COLS(c)
#endif

    public:

        CSeqInput_SAM( 
                const std::string & name, bool paired, 
                common::CFileBase::TCompression c );

        virtual int NCols( void ) const { return paired_ ? 2 : 1; }
        virtual bool Next( void );

        virtual const TData & Data( int col ) const
        {
            CHECK_COLS( col );
            if( col == 0 ) return d_0_;
            return d_1_;
        }

        virtual const TQual & Qual( int col ) const
        {
            CHECK_COLS( col );
            return q_[col];
        }

    private:

        static const unsigned int FLAG_SINGLE = 5;
        static const unsigned int FLAG_PAIRED_1 = 77;
        static const unsigned int FLAG_PAIRED_2 = 141;

        static const unsigned int ID_FIELD   = 0;
        static const unsigned int FLAG_FIELD = 1;
        static const unsigned int SEQ_FIELD  = 9;
        static const unsigned int QUAL_FIELD = 10;

        static const unsigned int MIN_FIELDS = 11;

        struct SEntry {
            std::string id;
            std::string seq;
            std::string qual;
            unsigned int flag;
            unsigned int idx;
        };

        SEntry ParseEntry( void );

        bool NextUnpaired( void );
        bool NextPaired( void );

        void ReadLine( void );
        void SkipHeader( void );
        void MakeWords( void );

        std::string CheckIDs( 
                const std::string & id_1, const std::string & id_2 );

        bool paired_;
        std::string name_;
        // std::auto_ptr< common::CReadTextFile > in_;
        std::unique_ptr< common::CReadTextFile > in_;

        TData d_0_, d_1_;
        TSeq s_[MAX_COLS];
        TSeqSize sz_[MAX_COLS];
        std::string q_[MAX_COLS];

        std::string line_;
        std::vector< std::string > words_;
};

END_NS( seq )
END_STD_SCOPES

#endif

