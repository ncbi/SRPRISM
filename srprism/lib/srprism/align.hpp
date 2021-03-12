/*  $Id: align.hpp 536631 2017-05-22 12:51:58Z morgulis $
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
 * File Description: srprism extension aligner
 *
 */

#ifndef __SRPRISM_ALIGN_HPP__
#define __SRPRISM_ALIGN_HPP__

#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif

#ifndef NCBI_CPP_TK

#include "../common/def.h"

#include <cstring>

#include "../common/bits.hpp"
#include "../seq/seqdef.hpp"
#include "srprismdef.hpp"
#include "stat.hpp"
#include "query_data.hpp"
#include "seqstore.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>

#include <cstring>

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/stat.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/query_data.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CBreakSegs
{
    private:

        typedef std::vector< TSeqSize > TVertices;
        typedef std::vector< int > TSiblings;
        typedef std::vector< common::Uint1 > TCounts;
        typedef std::vector< bool > TMarks;

    public:

        CBreakSegs()
        {
            left_.reserve( 2*MAX_N_HASHES );
            right_.reserve( 2*MAX_N_HASHES );
            siblings_.reserve( 2*MAX_N_HASHES );
            counts_.reserve( 2*MAX_N_HASHES );
            target_.reserve( 2*MAX_N_HASHES );
            marks_.reserve( 2*MAX_N_HASHES );
        }

        int Add( 
                TSeqSize start, TSeqSize end, 
                common::Uint1 target, int sibling = -1 )
        {
            left_.push_back( start );
            right_.push_back( end );
            siblings_.push_back( sibling );
            counts_.push_back( 0 );
            target_.push_back( target );
            marks_.push_back( false );
            return Size() - 1;
        }

        size_t Size( void ) const { return left_.size(); }

        bool Check( void ) const
        {
            for( int i( 0 ); i < (int)Size(); ++i ) {
                if( marks_[i] ) {
                    if( counts_[i] < target_[i] ) {
                        int j( siblings_[i] );

                        if( j != -1 ) {
                            if( !marks_[j] ) continue;
                            if( counts_[j] < target_[j]) return false;
                        }
                        else return false;
                    }
                }
            }

            return true;
        }

        void SetCurr( TSeqSize curr, bool d = false, bool l = true )
        {
            for( size_t i( 0 ); i < Size(); ++i ) {
                if( marks_[i] && left_[i] <= curr && right_[i] > curr ) {
                    if( d && ((left_[i] == curr && !l) || 
                              (right_[i] == curr + 1 && l)) ) {
                        if( counts_[i] == 0 ) {
                            if( l ) SetCurr( curr + 1 );
                            else SetCurr( curr - 1 );

                            continue;
                        }
                    }

                    ++counts_[i];
                }
            }
        }

        void SetBounds( TSeqSize start, TSeqSize end )
        {
            for( int i( 0 ); i < (int)Size(); ++i ) {
                marks_[i] = false; 
                counts_[i] = 0;
                if( left_[i] >= start && right_[i] <= end ) marks_[i] = true;
            }
        }

        void Clear( void )
        {
            left_.clear();
            right_.clear();
            siblings_.clear();
            counts_.clear();
            target_.clear();
            marks_.clear();
        }

        TSeqSize Left( int idx ) const { return left_[idx]; }
        TSeqSize Right( int idx ) const { return right_[idx]; }
        bool Empty() const { return (Size() == 0); }

    private:

        TVertices left_, right_;
        TSiblings siblings_;
        TCounts counts_, target_;
        TMarks marks_;
};

//------------------------------------------------------------------------------
namespace {
    template< typename int_t > struct CFAC_MASK
    { static const int_t VALUE = (int_t)0xaaaaaaaaaaaaaaaaULL; };

    template< typename int_t > struct LL_MASK
    { static const int_t VALUE = (int_t)0xfffffffffffffffcULL; };
}

template< typename int_t, int n_err >
struct CFastAlignCheck
{ bool operator()( int_t l, int_t r ) const { return false; } };

template< typename int_t, int n_err >
struct CFastAlignCheckWithMasks {};

template< typename int_t >
struct CFastAlignCheck< int_t, 1 >
{
    static const int_t MASK = CFAC_MASK< int_t >::VALUE;
    static const int_t LL_MASK = LL_MASK< int_t >::VALUE;

    bool operator()( int_t l, int_t r, int_t mask = LL_MASK ) const
    {
        if( l == r ) return false;

        // set even bits according to differing letters
        //
        int_t conv( l^r );
        conv |= (conv<<1);
        conv &= MASK;

        // remove the first set bit; 0 means one mismatch
        //
        int_t t( conv&(conv - 1) );
        if( t == 0 ) return true;

        // check for indel
        int b( common::FirstSetBit_Left< int_t >( conv ) );
        l <<= b;
        r <<= b;
        mask <<= b;
        if( (l<<2) == (r&mask) || (r<<2) == (l&mask) ) return true;

        return false;
    }
};

template< typename int_t >
struct CFastAlignCheckWithMasks< int_t, 1 >
{
    static const int_t MASK = CFAC_MASK< int_t >::VALUE;
    static const int_t LL_MASK = LL_MASK< int_t >::VALUE;

    bool operator()( 
            int_t l, int_t r, int_t m1, int_t m2, 
            int_t mask = LL_MASK ) const
    {
        // set even bits according to differing letters
        //
        int_t conv( l^r );
        conv |= (conv<<1);
        conv &= MASK;
        conv |= (m1|m2);

        // remove the first set bit; 0 means one mismatch
        //
        int_t t( conv&(conv - 1) );
        if( t == 0 ) return true;

        // check for indel
        int b( common::FirstSetBit_Left< int_t >( conv ) );
        l <<= b;
        r <<= b;
        m1 <<= b;
        m2 <<= b;
        mask <<= b;
        if( m2 == 0 && (m1<<2) == 0 && (l<<2) == (r&mask) ) return true;
        if( m1 == 0 && (m2<<2) == 0 && (r<<2) == (l&mask) ) return true;

        return false;
    }
};

template< typename int_t >
struct CFastAlignCheck< int_t, 2 >
{
    static const int_t MASK = CFAC_MASK< int_t >::VALUE;
    static const int_t LL_MASK = LL_MASK< int_t >::VALUE;

    bool operator()( int_t l, int_t r, int_t mask = LL_MASK ) const
    {
        if( l == r ) return false;
        if( CFastAlignCheck< int_t, 1 >()( l, r, mask ) ) return false;

        // set even bits according to differing letters
        //
        int_t conv( l^r );
        conv |= (conv<<1);
        conv &= MASK;

        // remove the first 2 set bits; 0 means two mismatch
        //
        int_t t( conv&(conv - 1) );
        t &= (t-1);
        if( t == 0 ) return true;

        // remove first error and check for the remainder for 
        // the possibility of 1-error alignment
        //
        int b( common::FirstSetBit_Left< int_t >( conv ) );
        l <<= b;
        r <<= b;
        mask <<= b;
        int_t l1( l&mask ), r1( r&mask );
        l <<= 2;
        r <<= 2;
        mask <<= 2;

        if( CFastAlignCheck< int_t, 1 >()( l, r1, mask ) ||
                CFastAlignCheck< int_t, 1 >()( l1, r, mask ) ||
                CFastAlignCheck< int_t, 1 >()( l, r, mask ) ) {
            return true;
        }

        return false;
    }
};

template< typename int_t >
struct CFastAlignCheckWithMasks< int_t, 2 >
{
    static const int_t MASK = CFAC_MASK< int_t >::VALUE;
    static const int_t LL_MASK = LL_MASK< int_t >::VALUE;

    bool operator()( 
            int_t l, int_t r, int_t m1, int_t m2, 
            int_t mask = LL_MASK ) const
    {
        // set even bits according to differing letters
        //
        int_t conv( l^r );
        conv |= (conv<<1);
        conv &= MASK;
        conv |= (m1|m2);

        // remove the first 2 set bits; 0 means two mismatches
        //
        int_t t( conv&(conv - 1) );
        t &= (t-1);
        if( t == 0 ) return true;

        // remove first error and check for the remainder for 
        // the possibility of 1-error alignment
        //
        int b( common::FirstSetBit_Left< int_t >( conv ) );
        l <<= b;
        r <<= b;
        m1 <<= b;
        m2 <<= b;
        mask <<= b;
        int_t l1( l&mask ), r1( r&mask );
        l <<= 2;
        r <<= 2;
        mask <<= 2;

        if( CFastAlignCheckWithMasks< int_t, 1 >()( 
                    l, r1, (m1<<2), m2, mask ) ||
                CFastAlignCheckWithMasks< int_t, 1 >()( 
                    l1, r, m1, (m2<<2), mask ) ||
                CFastAlignCheckWithMasks< int_t, 1 >()( 
                    l, r, (m1<<2), (m2<<2), mask ) ) {
            return true;
        }

        return false;
    }
};

//------------------------------------------------------------------------------
struct SMatrixEntry
{
    common::Uint4 serial;
    common::Uint2 range;
    common::Uint2 penalty;
};

//------------------------------------------------------------------------------
struct CExtensionSpaceAllocator
{
    typedef std::pair< SMatrixEntry *, common::Uint4 * > TExtensionSpaceHandle;

    CExtensionSpaceAllocator(
            common::Uint1 max_n_id,
            TSeqSize max_query_len,
            size_t max_hits )
        : serial_( 0xFFFFFFFF )
    {
        TSeqSize col_shift( common::BinLog( 2*max_n_id + 3 ) + 1 ),
                 col_len( 1<<col_shift );
        size_t pool_sz( col_len*max_query_len*max_hits );
        m_pool_.resize( pool_sz );
        q_pool_.resize( pool_sz );
    }

    void Clean( void )
    {
        if( ++serial_ == 0 ) {
            memset( &m_pool_[0], 0, sizeof( SMatrixEntry )*m_pool_.size() );
            ++serial_;
        }
    }

    common::Uint4 GetSerial( void ) const { return serial_; }

    TExtensionSpaceHandle Alloc( TSeqSize qlen, common::Uint1 n_diag )
    {
        SRPRISM_ASSERT( n_diag*qlen < m_pool_.size() );
        TExtensionSpaceHandle res( std::make_pair( &m_pool_[0], &q_pool_[0] ) );
        return res;
    }

    private:

        common::Uint4 serial_;
        std::vector< SMatrixEntry > m_pool_;
        std::vector< common::Uint4 > q_pool_;
};

//------------------------------------------------------------------------------
class CHit
{
    private:

        typedef common::Uint8 TDWord;

        static const size_t LBITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;
        static const size_t WLETTERS = 
            sizeof( TWord )*seq::SCodingTraits< SEQDATA_CODING >::PACK_FACTOR;
        static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;
        static const size_t WMASK = 
            common::SBitFieldTraits< size_t, WSHIFT >::MASK;
        static const size_t WBITS = WLETTERS*LBITS;
        static const TDWord LWORD_MASK = 
            common::SBitFieldTraits< TDWord, WBITS >::MASK;

        static const common::Uint2 ERR_STEP = 0x100;

        typedef CSeqStore::TPos TPos;

    public:

        struct SErrRec
        {
            common::Uint2 qpos;
            TErrType err_type;
        };

    private:

        struct STraceBack
        {
            STraceBack() : mat( 0 ) {}
            bool IsEmpty( void ) const { return (mat == 0); }
            int NErr( void ) const { return (int)(mat[best_idx].penalty>>8); }
            int NId( void ) const { return (int)(mat[best_idx].penalty&0xFF); }

            void Reset( void ) { start_idx = best_idx; }

            bool NextError( SErrRec & result )
            {
                const SMatrixEntry * e( mat + start_idx );
                if( e->penalty == 0 ) return false;
                result.qpos = start_idx/col_len;
                const SMatrixEntry * c;

                if( e - col_len >= mat ) {
                    c = e - col_len;

                    if( c->serial == serial ) {
                        c -= (c->range<<col_shift);
    
                        if( c->serial == serial && 
                                c->penalty + ERR_STEP == e->penalty ) {
                            result.err_type = SErrType::M;
                            start_idx = c - mat;
                            return true;
                        }
                    }
                }

                c = e + 1;

                if( c->serial == serial ) {
                    c -= (c->range<<col_shift);

                    if( c->serial == serial && 
                            c->penalty + ERR_STEP + 1 == e->penalty ) {
                        result.err_type = SErrType::D;
                        start_idx = c - mat;
                        return true;
                    }
                }

                if( e - col_len - 1 >= mat ) {
                    c = e - col_len - 1;

                    if( c->serial == serial ) {
                        c -= (c->range<<col_shift);

                        if( c->serial == serial && 
                                c->penalty + ERR_STEP + 1 == e->penalty ) {
                            result.err_type = SErrType::I;
                            start_idx = c - mat;
                            return true;
                        }
                    }
                }

                // control flow should never get here
                //
                SRPRISM_ASSERT( false );
                return false; 
            }

            template< int dir, bool partial_align > void ExtendNonExact(
                const TWord * q, const TWord * m,
                const TWord * s, const TWord * sm,
                TSeqSize qwoff, TSeqSize swoff,
                TSeqSize q_min, TSeqSize q_max, TSeqSize q_lim, TSeqSize s_lim, 
                common::Uint1 n_err, common::Uint1 n_id,
                CExtensionSpaceAllocator & ma );

            void TrimAndUpdate( 
                    common::Uint4 curr_idx, 
                    TSeqSize q_lim, TSeqSize q_min, TSeqSize q_max,
                    TSeqSize range, int n_err );

            common::Uint4 start_idx, best_idx, serial;
            int best_score;
            TSeqSize col_len, col_shift, flank;
            SMatrixEntry * mat;
        };

        template< int dir > struct SExtendDirTraits {};

        static inline TWord CutWord( TDWord dword, const TSeqSize & offset );

        static inline void ExtendLeftExact(
                const TWord * q, const TWord * qm,
                const TWord * s, const TWord * sm,
                TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
                TDWord dq, TDWord dm, TDWord ds, TDWord dsm,
                TSeqSize & n_matched );

        static inline void ExtendRightExact(
                const TWord * q, const TWord * qm,
                const TWord * s, const TWord * sm,
                TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
                TDWord dq, TDWord dm, TDWord ds, TDWord dsm,
                TSeqSize & n_matched );

    public:

        typedef std::pair< common::Uint2, common::Uint2 > TBreakSeg;
        typedef std::vector< TBreakSeg > TBreakSegs;
        typedef std::vector< SErrRec > TErrors;

        CHit() : qdata_( 0, 0, 0 ) {}

        CHit( const CQueryData & q, TStrand s )
            : qdata_( q ), anchor_( 0 ), left_flank_( 0 ), align_len_( 0 ), 
              strand_( s ), n_id_( qdata_.NIns() + qdata_.NDel() ),
              n_d_( qdata_.NDel() )
        {}

        template< bool partial_alignment >
        bool Extend( 
                const CSeqStore & ss, CExtensionSpaceAllocator & ma,
                TSeqSize sa_start, TSeqSize sa_end,
                TSeqSize seed_qoff, TPos seed_soff,
                TSeqSize seed_slen, int n_err,
                CBreakSegs * bsegs, size_t * stat, size_t * ustat );

        CQueryData GetQueryData( void ) const { return qdata_; }
        int FullNId( void ) const { return n_id_; } 
        int FullNDel( void ) const { return n_d_; }
        int FullNErr( void ) const { return err_data_.size(); }
        int FullNIns( void ) const { return n_id_ - n_d_; }
        int FullNGOpen( void ) const { return 0; }
        TStrand Strand( void ) const { return strand_; }
        TSeqSize AlignLen( void ) const { return align_len_; }
        TSeqSize LeftOffset( void ) const { return left_flank_; }

        TSeqSize RightOffset( void ) const 
        { return qdata_.Len() - left_flank_ - AlignLen() + n_d_; }

        TSeqSize AnchorSubjLen( void ) const
        { return AlignLen() + n_d_ - n_id_; }

        TSeqSize SubjLen( void ) const 
        { return AlignLen() + n_d_ - n_id_ + LeftOffset() + RightOffset(); }

        TPos Anchor( void ) const { return anchor_; }

        const TErrors & GetErrors( void ) const { return err_data_; }

    private:

        CQueryData qdata_;
        TPos anchor_;
        TSeqSize left_flank_, align_len_;
        TStrand strand_;
        TErrors err_data_;
        int n_id_, n_d_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

