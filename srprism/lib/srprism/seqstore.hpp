/*  $Id: seqstore.hpp 431273 2014-04-02 17:10:44Z morgulis $
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
 * File Description: part of the index containing packed subject sequence data
 *
 */

#ifndef __SRPRISM_SEQSTORE_HPP__
#define __SRPRISM_SEQSTORE_HPP__

#include "../common/def.h"

#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <common/bits.hpp>
#include <seq/seqdef.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/seqstore_base.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore_base.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CSeqStore : public CSeqStoreBase
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode VALID = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == VALID ) return "sequence store validation error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CSeqStore( const std::string & basename, CMemoryManager & mem_mgr );
        ~CSeqStore() { Unload(); }

        void Load(void);
        void Unload( void );
        void LoadAmbigData( void );
        void UnloadAmbigData( void );
        common::Uint4 OverlapFactor( void ) const { return max_seq_overlap_; }

    private:

        CSeqStore( const CSeqStore & );
        CSeqStore & operator=( const CSeqStore & );

        typedef std::vector< TMaskUnit > TAmbigMask;

        struct SSeqMapEntry
        {
            TDBOrdId oid;
            TPos seq_start;
            TPos body_start;
            TPos body_end;
            TPos seq_end;
            TSeqSize ref_loc_start;
            TSeqSize ref_loc_end;

            friend bool operator<( 
                    const SSeqMapEntry & l, const SSeqMapEntry & r )
            { return l.seq_start < r.seq_start; }

            friend std::ostream & operator<<( 
                    std::ostream & os, SSeqMapEntry const & e ) {
                os << "oid: " << e.oid << "; "
                   << "sstart: " << e.seq_start << "; "
                   << "bstart: " << e.body_start << "; "
                   << "bend: " << e.body_end << "; "
                   << "send: " << e.seq_end << "; "
                   << "rstart: " << e.ref_loc_start << "; "
                   << "rend: " << e.ref_loc_end;
                return os;
            }
        };

        typedef std::vector< SSeqMapEntry > TSeqMap;
        typedef std::vector< TDBOrdId > TALLists;

        struct SRefSegEnd
        {
            size_t all_idx;
            TDBOrdId ref_sid;
            TSeqSize ref_pos;

            SRefSegEnd( size_t a, TDBOrdId s, TSeqSize p )
                : all_idx( a ), ref_sid( s ), ref_pos( p )
            {}

            friend bool operator<( const SRefSegEnd & l, const SRefSegEnd & r )
            {
                if( l.ref_sid == r.ref_sid ) return l.ref_pos < r.ref_pos;
                return l.ref_sid < r.ref_sid;
            }
        };

        typedef std::vector< SRefSegEnd > TRefSegs;

        void LoadHeader( void );
        void LoadSeqMap( void );
        void ComputeSeqOverlap( void );
        void LoadDynamicData( void );

        std::string basename_;
        common::Uint4 n_seq_;
        common::Uint4 data_sz_;
        common::Uint4 ambig_map_sz_;
        common::Uint4 ambig_data_sz_;
        TAmbigMask ambig_mask_;
        TSeqMap seq_map_;
        TALLists al_lists_;
        TRefSegs ref_segs_;
        CMemoryManager & mem_mgr_;
        SAmbigRun * ambig_map_;
        TLetter * ambig_data_;
        TWord * seq_data_;
        TWord * rev_seq_data_;
        size_t segment_letters_;

        common::Uint4 max_seq_overlap_;

    public:

        typedef std::pair< TDBOrdId, TSeqSize > TDecSeqData;

        TDecSeqData DecodePos( TPos pos ) const
        {
            SSeqMapEntry t; t.seq_start = pos;
            TSeqMap::const_iterator p( std::upper_bound( 
                        seq_map_.begin(), seq_map_.end(), t ) );
            SRPRISM_ASSERT( p != seq_map_.begin() );
            --p;
            return std::make_pair( p - seq_map_.begin(), pos - p->seq_start );
        }

        TPos EncodePos( TDecSeqData pos ) const
        { return seq_map_[pos.first].seq_start + pos.second; }

        common::Sint8 GetAdjustedPos( TDBOrdId sid, TSeqSize pos ) const
        {
            const SSeqMapEntry & s( seq_map_[sid] );
            return (common::Sint8)pos - (s.body_start - s.seq_start);
        }

        size_t NSeq( void ) const { return n_seq_; }

        std::pair< const TWord *, TSeqSize > FwDataPtr( TPos pos ) const
        {
            return std::make_pair( 
                    seq_data_ + (pos>>WORD_SHIFT), (pos&WORD_MASK) );
        }

        std::pair< const TWord *, TSeqSize > RvDataPtr( TPos pos ) const
        {
            pos = data_sz_ - pos;
            return std::make_pair(
                    rev_seq_data_ + (pos>>WORD_SHIFT), (pos&WORD_MASK) );
        }

        TDBOrdId GetRefOId( TDBOrdId oid ) const
        {
            SRPRISM_ASSERT( oid < seq_map_.size() );
            return seq_map_[oid].oid;
        }

        TSeqSize FwTailLen( TPos pos ) const
        {
            SSeqMapEntry t; t.seq_start = pos;
            TSeqMap::const_iterator it( 
                    std::upper_bound( seq_map_.begin(), seq_map_.end(), t ) );
            SRPRISM_ASSERT( it != seq_map_.begin() );
            --it;
            return it->seq_end - pos;
        }

        TSeqSize RvTailLen( TPos pos ) const
        {
            SSeqMapEntry t; t.seq_start = pos;
            TSeqMap::const_iterator it( 
                    std::upper_bound( seq_map_.begin(), seq_map_.end(), t ) );
            SRPRISM_ASSERT( it != seq_map_.begin() );
            --it;
            return pos - it->seq_start;
        }

        TSeqSize GetSeqLen( TDBOrdId s_n ) const
        { 
            SRPRISM_ASSERT( s_n < seq_map_.size() );
            return seq_map_[s_n].seq_end - seq_map_[s_n].seq_start; 
        }

        TPos NextAmbigPos( TPos pos ) const
        {
            if( ambig_map_sz_ == 0 ) return (TPos)data_sz_;
            SAmbigRun target( pos, 0 );
            SAmbigRun * loc( std::upper_bound( 
                        ambig_map_, ambig_map_ + ambig_map_sz_, target ) ),
                      * loc_1( loc - 1 );
            
            if( loc != ambig_map_ ) {
                if( loc_1->pos + loc_1->len <= pos ) {
                    if( loc != ambig_map_ + ambig_map_sz_ ) return loc->pos;
                    else return (TPos)data_sz_;
                }
                else return pos;
            }
            else return loc->pos;
        }

        TLetter GetLetter( TDBOrdId snum, TSeqSize spos ) const
        {
            TPos pos( EncodePos( std::make_pair( snum, spos ) ) );
            if( NextAmbigPos( pos ) == pos ) return (TLetter)'N';
            std::pair< const TWord *, TSeqSize > wd( FwDataPtr( pos ) );
            TLetter l( seq::GetStreamLetter< SEQDATA_CODING >( 
                        wd.first, wd.second ) );
            return seq::SCodingTraits< 
                SEQDATA_CODING >::Recode< OUTPUT_CODING >( l );
        }

        void GetInsideList( 
                TDBOrdId ref_oid, TSeqSize start,
                std::vector< TDBOrdId > & oid_list ) const
        {
            SRefSegEnd t( 0, ref_oid, start );
            TRefSegs::const_iterator i( std::upper_bound( 
                        ref_segs_.begin(), ref_segs_.end(), t ) );

            if( i != ref_segs_.begin() ) {
                size_t n;
                
                if( i == ref_segs_.end() ) {
                    n = al_lists_.size() - (--i)->all_idx;
                }
                else {
                    --i;
                    n = (i+1)->all_idx - i->all_idx;
                }

                for( size_t j( 0 ); j < n; ++j ) {
                    oid_list.push_back( al_lists_[i->all_idx + j] );
                }
            }
        }

    private:

        size_t NAmbigs_Slow( TPos pos1, TPos pos2 ) const
        {
            size_t res( 0 );

            for( TPos pos( pos1 ); pos < pos2; ) {
                TPos pos_1( NextAmbigPos( pos ) );

                if( pos_1 == pos ) { ++res; ++pos; }
                else pos = pos_1;
            }

            return res;
        }

        bool SegMask( size_t seg ) const
        {
            size_t unit( seg>>MASK_SHIFT ), bit( seg&MASK_MASK );
            return common::GetBit( bit, ambig_mask_[unit] );
        }

    public:

        size_t NAmbigs( TPos pos1, TPos pos2 ) const
        {
            size_t res( 0 );
            size_t segment_shift( SegmentShift( segment_letters_ ) );

            for( TPos pos( pos1 ); pos < pos2; ) {
                size_t seg( pos>>segment_shift );
                TPos pos_e( std::min( pos2, (TPos)((seg+1)<<segment_shift) ) );
                if( SegMask( seg ) ) res += NAmbigs_Slow( pos, pos_e );
                pos = pos_e;
            }

            return res;
        }

        bool CheckRegion( TPos start_pos, TSeqSize len ) const
        {
            TDBOrdId oid( DecodePos( start_pos ).first );
            const SSeqMapEntry & s( seq_map_[oid] );
            return s.body_start < start_pos + len && start_pos < s.body_end;
        }

        bool CheckRegionPair( 
                TPos start_pos_1, TPos start_pos_2, 
                TSeqSize len_1, TSeqSize len_2 ) const
        { 
            return CheckRegion( start_pos_1, len_1 ) || 
                   CheckRegion( start_pos_2, len_2 );
        }
};

END_NS( srprism )
END_STD_SCOPES

#endif

