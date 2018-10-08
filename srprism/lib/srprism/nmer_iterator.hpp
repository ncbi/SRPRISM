/*  $Id: nmer_iterator.hpp 205515 2010-09-20 15:51:04Z morgulis $
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
 * File Description: sequence traversal by conseequtive n-mers
 *
 */

#ifndef __PDOP_NMER_ITERATOR_HPP__
#define __PDOP_NMER_ITERATOR_HPP__

#include "../common/def.h"

#include "../seq/seqdef.hpp"
#include "srprismdef.hpp"
#include "seqstore.hpp"

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CNMerIterator
{
    public:

        CNMerIterator( const CSeqStore & ss, TDBOrdId seq_n );
        
        bool End( void ) const { return end_; }
        TSeqSize Pos( void ) const { return curr_pos_ + HASH_LEN - 1; }

        TPrefix Prefix( void ) const 
        { 
            return curr_prefix_ <= reverse_prefix_ ? curr_prefix_ 
                                                   : reverse_prefix_; 
        }

        TStrand Strand( void ) const
        { 
            return curr_prefix_ <= reverse_prefix_ ? seq::STRAND_FW 
                                                   : seq::STRAND_RV; 
        }

        bool Palindrome( void ) const 
        { return curr_prefix_ == reverse_prefix_; }

        bool Next( void )
        { 
            if( !End() ) NextLetter();
            while( !End() && len_ < HASH_LEN ) NextLetter(); 
            return !End();
        }

    private:

        bool NextLetter( void );
        const CSeqStore & ss_;
        TSeqSize len_, left_, fw_off_;
        CSeqStore::TPos curr_pos_;
        const TWord * fw_ptr_;
        TPrefix curr_prefix_, reverse_prefix_;
        bool end_;
};


END_NS( srprism )
END_STD_SCOPES

#endif

