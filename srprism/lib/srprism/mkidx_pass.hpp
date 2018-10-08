/*  $Id: mkidx_pass.hpp 210557 2010-11-04 17:15:08Z morgulis $
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
 * File Description: one pass of index creation
 *
 */

#ifndef __SRPRISM_MKIDX_PASS_HPP__
#define __SRPRISM_MKIDX_PASS_HPP__

#include "../common/def.h"

#include "../common/binfile.hpp"
#include "../common/bits.hpp"
#include "../common/trace.hpp"
#include "seqstore.hpp"
#include "index_base.hpp"

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CMkIdxPass : public CIndexBase
{
    private:

        static const size_t FORMAT_THRESHOLD = 100;

        typedef CSeqStore::TPos TPos;

        struct SPosEntry
        {
            TPos pos;
            TStrand strand;

            union {
                TPrefix prefix;
                TExtension extension;
            };

            friend bool operator<( 
                    const SPosEntry & lhs,
                    const SPosEntry & rhs )
            {
                return (lhs.prefix == rhs.prefix) ? 
                       (lhs.strand == rhs.strand) ? (lhs.pos < rhs.pos)
                                                  : (lhs.strand < rhs.strand)
                                                  : (lhs.prefix < rhs.prefix);
            }
        };

        class CMapPrefixIterator
        {
            public:

                static const size_t MAP_PREFIX_SIZE = (size_t)MAP_PREFIX_LEN;
                static const size_t MASK_BITS =
                    seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS*(
                            PREFIX_LEN - MAP_PREFIX_SIZE);
                static const size_t MASK = 
                    common::SBitFieldTraits< size_t, MASK_BITS >::MASK;

                CMapPrefixIterator( 
                        const SPosEntry * entries, size_t n_entries )
                    : entries_( entries ), n_entries_( n_entries ),
                      start_entry_( 0 ), end_entry_( 0 ), 
                      last_nmer_( entries_[n_entries_ - 1].prefix )
                { 
                    assert( n_entries_ > 0 );
                    Next(); 
                }

                bool Next( void )
                {
                    if( End() ) return false;

                    for( start_entry_ = end_entry_; 
                            end_entry_ < n_entries_ &&
                            StartPrefix() == EndPrefix(); ++end_entry_ );

                    return true;
                }

                bool End( void ) const { return start_entry_ == n_entries_; }

                size_t StartEntry( void ) const { return start_entry_; }
                size_t EndEntry( void ) const { return end_entry_; }

                size_t StartPrefix( void ) const 
                { return (entries_[start_entry_].prefix>>MASK_BITS); }

                size_t EndPrefix( void ) const 
                { 
                    if( end_entry_ != n_entries_ ) {
                        return (entries_[end_entry_].prefix>>MASK_BITS); 
                    }
                    else return ((last_nmer_ + 1)>>MASK_BITS);
                }

            private:

                const SPosEntry * entries_;
                size_t n_entries_;
                size_t start_entry_, end_entry_;
                size_t last_nmer_;
        };

        class CPrefixIterator
        {
            public:

                CPrefixIterator( 
                        const SPosEntry * entries,
                        size_t start, size_t end )
                    : entries_( entries ), start_( start ), end_( end ),
                      cend_( start )
                { Next(); }

                bool Next( void )
                {
                    if( End() ) return false;

                    for( start_ = cend_; 
                            cend_ != end_ && 
                            entries_[cend_].prefix == entries_[start_].prefix;
                            ++cend_ );

                    return true;
                }

                size_t StartEntry( void ) const { return start_; }
                size_t EndEntry( void ) const { return cend_; }

                bool End( void ) const { return start_ == end_; }

            private:

                const SPosEntry * entries_;
                size_t start_, end_, cend_;
        };

        typedef CPrefixIterator TExtensionIterator;

    public:

        static const TSeqSize HASH_KEY_SIZE = 8;

    private:

        static const size_t LETTER_BITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t HASH_KEY_BITS = HASH_KEY_SIZE*LETTER_BITS;
        static const size_t NUM_HASH_KEYS = 
            common::SBitFieldTraits< size_t, HASH_KEY_BITS >::MAX + 1;
        static const size_t SHIFT = 
            common::BYTEBITS*sizeof( TPrefix ) - HASH_KEY_BITS;

    public:


        CMkIdxPass( const size_t * counts_table,
                    const CSeqStore & seq_store, 
                    void * free_space, size_t free_space_size, size_t & start_idx, 
                    common::CWriteBinFile & map_file,
                    common::CWriteBinFile & rmap_file,
                    common::CWriteBinFile & idx_file );

        ~CMkIdxPass() 
        { 
            M_TRACE( common::CTracer::INFO_LVL, 
                     "repeat map size is " << rmap_size_ ); 
        }

        void Run( void );

    private:

        static std::pair< size_t, size_t > SplitByStrand(
                SPosEntry * entries, size_t start, size_t end );

        size_t Estimate( size_t start_entry, size_t end_entry );
        size_t EstimateNormal( size_t start, size_t end );
        size_t EstimateSpecial( size_t start, size_t end );

        void UpdateRMap( size_t start, size_t end );

        size_t CountExtensions( 
                SPosEntry * entries, size_t start, size_t end ) const;
        size_t EstimateExtensions( 
                SPosEntry * entries, size_t start, size_t end ) const;

        void SetUpPalindromeData( 
                std::vector< SPosEntry > & entries, size_t start, size_t end );
        void SetUpExtensions( 
                SPosEntry * entries, size_t start, size_t end, TStrand s );
        void FlushMapPrefix( size_t start_entry, size_t end_entry );
        void FlushSpecialForStrand( 
                SPosEntry * entries, size_t start, size_t end );
        void FlushSpecial( size_t start, size_t end );
        void FlushNormal( size_t start, size_t end, common::Uint2 descr );

        const CSeqStore & seq_store_;
        common::CWriteBinFile & map_file_;
        common::CWriteBinFile & rmap_file_;
        common::CWriteBinFile & idx_file_;
        size_t start_idx_, end_idx_;
        size_t total_entries_;
        SPosEntry * entries_;
        size_t rmap_size_;
};

END_NS( srprism )
END_STD_SCOPES

#endif
