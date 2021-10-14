/*  $Id: bnf.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: Bad N-mer filter.
 *
 */

#ifndef __SRPRISM_BNF_HPP__
#define __SRPRISM_BNF_HPP__

#ifndef NCBI_CPP_TK

#include "../common/def.h"
#include "srprismdef.hpp"
#include "align.hpp"

#else

#include <../src/internal/align_toolbox/srprism/lib/common/def.h>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/align.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
namespace {
    template< typename t_uint >
    static t_uint RemoveLetter_R( t_uint v, size_t idx )
    {
        static const size_t LBITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t LSHIFT = SBinLog< LBITS >::VALUE;
        static const size_t NBITS = sizeof( t_uint )*BYTEBITS;
        static const size_t NLETTERS = (NBITS>>LSHIFT);

        idx = ((NLETTERS - idx - 1)<<LSHIFT);
        t_uint mask( MaskBack< t_uint >( idx ) );
        return ((v>>LBITS)&~mask) + (v&mask);
    }

    template< typename t_uint, int n_err >
    struct SRemoveLettersIterator
    {
        static const size_t LBITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t LSHIFT = SBinLog< LBITS >::VALUE;
        static const size_t NBITS = sizeof( t_uint )*BYTEBITS;
        static const size_t NLETTERS = (NBITS>>LSHIFT);

        SRemoveLettersIterator() 
            : word_( 0 ), reduced_word_( 0 ), letter_idx_( NLETTERS ) {}

        SRemoveLettersIterator( t_uint word ) 
            : word_( word ), 
              reduced_word_( RemoveLetter_R< t_uint >( word, 0 ) ),
              letter_idx_( 0 )
        {
            rec_iter_.Init( reduced_word_, 1 );
        }

        void Init( t_uint word, size_t idx = 0 )
        { 
            word_ = word; 
            reduced_word_ = RemoveLetter_R< t_uint >( word, 0 );
            letter_idx_ = idx; 
        }

        bool Next( t_uint & res )
        {
            if( letter_idx_ < NLETTERS - n_err + 1 ) {
                if( !rec_iter_.Next( res ) ) {
                    TLetter l( seq::GetLetter< SEQDATA_CODING, t_uint >(
                                word_, letter_idx_ ) );

                    while( letter_idx_ < NLETTERS - n_err &&
                                l == seq::GetLetter< SEQDATA_CODING, t_uint >( 
                                    word_, letter_idx_ + 1 ) ) {
                        ++letter_idx_;
                    }

                    if( letter_idx_ >= NLETTERS - n_err ) return false;

                    reduced_word_ = 
                        RemoveLetter_R< t_uint >( word_, ++letter_idx_ );
                    rec_iter_.Init( reduced_word_, letter_idx_ + 1 );
                    rec_iter_.Next( res );
                }

                return true;
            }
            else return false;
        }

        t_uint word_, reduced_word_;
        size_t letter_idx_;
        SRemoveLettersIterator< t_uint, n_err - 1 > rec_iter_;
    };

    template< typename t_uint >
    struct SRemoveLettersIterator< t_uint, 1 >
    {
        static const size_t LBITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
        static const size_t LSHIFT = SBinLog< LBITS >::VALUE;
        static const size_t NBITS = sizeof( t_uint )*BYTEBITS;
        static const size_t NLETTERS = (NBITS>>LSHIFT);

        SRemoveLettersIterator() 
            : word_( 0 ), reduced_word_( 0 ), letter_idx_( NLETTERS ) {}

        SRemoveLettersIterator( t_uint word )
            : word_( word ), 
              reduced_word_( RemoveLetter_R< t_uint >( word, 0 ) ), 
              letter_idx_( 0 )
        {}

        void Init( t_uint word, size_t idx = 0 ) 
        { 
            word_ = word; 
            reduced_word_ = RemoveLetter_R< t_uint >( word, idx );
            letter_idx_ = idx ; 
        }

        bool Next( t_uint & res )
        {
            if( letter_idx_ < NLETTERS ) {
                res = reduced_word_;
                TLetter l( seq::GetLetter< SEQDATA_CODING, t_uint >(
                            word_, letter_idx_ ) );

                while( letter_idx_ < NLETTERS - 1 && 
                            l == seq::GetLetter< SEQDATA_CODING, t_uint >( 
                                word_, letter_idx_ + 1 ) ) {
                    ++letter_idx_;
                }

                reduced_word_ = 
                    RemoveLetter_R< t_uint >( word_, ++letter_idx_ );
                return true;
            }
            else return false;
        }

        t_uint word_, reduced_word_;
        size_t letter_idx_;
    };
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
class CBadNMerFilter
{
    struct SQCompare
    {
        bool operator()( const t_qdata & lhs, const t_qdata & rhs )
        { return lhs.ext < rhs.ext; }
    };

    public:

        static const size_t MAX_QDATA_SIZE = 255;

        typedef std::vector< size_t > TQExtSet;

        typedef struct
        {
            ssize_t s_ext;
            TQExtSet q_ext_set;
        } TMatchSet;

        CBadNMerFilter( size_t & stat, bool randomize );

        void Init( 
                const t_sdata * sdata, const t_qdata * qdata,
                size_t sdata_sz, size_t qdata_sz,
                TSeqSize ext_len, int filter_n_err );

        // returns the min number of errors required to match
        // match set pairs; -1 on end
        //
        int Next( void );

        // returns data for a single subject extension
        //
        const TMatchSet & GetMatchSet( void ) const { return match_set_; }

        int GetNErr() const { return curr_n_err_; }

    private:

        static const int LBITS = 
            seq::SCodingTraits< SEQDATA_CODING >::LETTER_BITS;

        static const int MAX_FILTER_N_ERR = 2;

        // for each query extension we have a left and right 8-mer;
        // by removing up to 2 letters from an 8-mer we can get up
        // 28 different 6-mers; so the blowup factor is 56
        //
        static const size_t QEXT_INDICES_SIZE = 56*MAX_QDATA_SIZE;

        struct SState
        {
            static const int ALL     = 0;
            static const int EXACT   = 1;
            static const int ONE_ERR = 2;
            static const int TWO_ERR = 3;
            static const int DONE    = 4;
        };

        struct SLUTable
        {
            typedef std::vector< common::Uint1 > TSerial;
            typedef std::vector< common::Uint1 > TLen;
            typedef std::vector< common::Uint2 > TIdx;

            // look up table is indexed by N-mer with N = 7 for
            // 1 error filter; N = 6 for 2 error filter
            //
            static const size_t LUT_SIZE_1 = 16*KILOBYTE;
            static const size_t LUT_SIZE_2 = 4*KILOBYTE;

            TSerial serial;
            TLen len_l, len_r;
            TIdx idx_l, idx_r;
        };

        void CheckSerial()
        {
            bool ok( true );

            if( lut_.serial.size() != SLUTable::LUT_SIZE_1 ) ok = false;
            else
            {
                for( auto s : lut_.serial )
                {
                    if( s >= serial_ ) { ok = false; break; }
                }
            }

            if( !ok )
            {
                throw std::runtime_error( "BAD SERIAL" );
            }
        }

        void DumpLUT( void )
        {
            for( size_t i( 0 ); i < SLUTable::LUT_SIZE_1; ++i ) {
                if( lut_.serial[i] == serial_ ) {
                    std::cerr << "left  " << std::hex << (int)i << std::dec
                              << " idx: " << (int)lut_.idx_l[i]
                              << " len: " << (int)lut_.len_l[i];

                    for( size_t j( 0 ); j < lut_.len_l[i]; ++j ) {
                        std::cerr << ' ' << (int)(q_ext_indices_[lut_.idx_l[i] + j]);
                    }

                    std::cerr << std::endl;
                    std::cerr << "right " << std::hex << (int)i << std::dec
                              << " idx: " << (int)lut_.idx_r[i]
                              << " len: " << (int)lut_.len_r[i];

                    for( size_t j( 0 ); j < lut_.len_r[i]; ++j ) {
                        std::cerr << ' ' << (int)(q_ext_indices_[lut_.idx_r[i] + j]);
                    }

                    std::cerr << std::endl;
                }
            }
        }

        size_t GetLeftIdx( THalfWord hw ) const
        { return lut_.serial[hw] == serial_ ? lut_.idx_l[hw] : 0; }

        size_t GetRightIdx( THalfWord hw ) const
        { return lut_.serial[hw] == serial_ ? lut_.idx_r[hw] : 0; }

        size_t GetLeftLen( THalfWord hw ) const
        { return lut_.serial[hw] == serial_ ? lut_.len_l[hw] : 0; }

        size_t GetRightLen( THalfWord hw ) const
        { return lut_.serial[hw] == serial_ ? lut_.len_r[hw] : 0; }

        template< int n_err >
        size_t GetCombinedLenLeft( THalfWord hw ) const
        {
            size_t res( 0 );
            SRemoveLettersIterator< THalfWord, n_err > RLI( hw );
            THalfWord w( 0 );
            while( RLI.Next( w ) ) res += GetLeftLen( w );
            return res;
        }

        template< int n_err >
        size_t GetCombinedLenRight( THalfWord hw ) const
        {
            size_t res( 0 );
            SRemoveLettersIterator< THalfWord, n_err > RLI( hw );
            THalfWord w( 0 );
            while( RLI.Next( w ) ) res += GetRightLen( w );
            return res;
        }

        typedef common::Uint1 TQExtIdx;
        typedef std::vector< TQExtIdx > TQExtIndices;

        typedef std::vector< common::Uint4 > TDupChecks;

        int NextAll( void );
        int NextExact( void );
        int NextOneErr( void );
        int NextTwoErr( void );

        int NextExactRand( void );

        template< int n_err > void NextPriv( void );

        template< int n_err > 
        void ProcessLUTRow( TQExtIdx * start, Uint1 & len );

        template< int n_err > void ProcessLUTEntry( TQExtIdx q_idx );

        bool SetUpLUT( int state );
        template< int n_err > bool SetUpLUT( void );

        void SetState( int new_state )
        {
            s_ext_idx_ = -1;
            state_ = new_state;
            curr_n_err_ = 0;

            if( state_ == SState::ONE_ERR || state_ == SState::TWO_ERR ) {
                if( ++serial_ == 0 ) { // reinit relevant data
                    lut_.serial.assign( SLUTable::LUT_SIZE_1, 0 );
                    ++serial_;
                }

                dup_checks_.assign( qdata_sz_, 0xFFFFFFFFUL );
                curr_n_err_ = state_ - 1;
            }
        }

        const t_sdata * sdata_; // assumed sorted by extension value
        const t_qdata * qdata_; // assumed sorted by extension value
        size_t sdata_sz_;
        size_t qdata_sz_;
        TSeqSize ext_len_;
        int filter_n_err_,
            curr_n_err_;

        common::Uint1 serial_;
        int state_;

        TMatchSet match_set_;
        ssize_t s_ext_idx_;

        size_t curr_q_ext_;
        TWord s_extension_, ext_mask_;

        SLUTable lut_;
        TQExtIndices q_ext_indices_;
        TDupChecks dup_checks_;

        bool randomize_;
        
        size_t & stat_;

        std::vector< Uint4 > rnd_map_;
};

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
CBadNMerFilter< t_sdata, t_qdata >::CBadNMerFilter( 
        size_t & stat, bool randomize )
    : sdata_( 0 ), qdata_( 0 ), serial_( 0 ), randomize_( randomize ), 
      stat_( stat )
{
    match_set_.q_ext_set.reserve( MAX_QDATA_SIZE );

    lut_.serial.resize( SLUTable::LUT_SIZE_1, 0 );
    lut_.len_l.resize( SLUTable::LUT_SIZE_1, 0 );
    lut_.len_r.resize( SLUTable::LUT_SIZE_1, 0 );
    lut_.idx_l.resize( SLUTable::LUT_SIZE_1, 0 );
    lut_.idx_r.resize( SLUTable::LUT_SIZE_1, 0 );

    q_ext_indices_.resize( QEXT_INDICES_SIZE, 0 );
    dup_checks_.resize( MAX_QDATA_SIZE );
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
bool CBadNMerFilter< t_sdata, t_qdata >::SetUpLUT( int state )
{
    if( state == SState::ONE_ERR ) return SetUpLUT< 1 >();
    else if( state == SState::TWO_ERR ) return SetUpLUT< 2 >();
    else SRPRISM_ASSERT( false );

    return false;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
template< int n_err >
bool CBadNMerFilter< t_sdata, t_qdata >::SetUpLUT( void )
{
    static const size_t SHIFT = BYTEBITS*sizeof( THalfWord );
    static const TWord MASK = SBitFieldTraits< TWord, SHIFT >::MASK;

    typedef SRemoveLettersIterator< THalfWord, n_err > TRLI;

    {
        size_t cnt( 0 );

        for( size_t i( 0 ); i < qdata_sz_; ++i ) {
            if( qdata_[i].count > 0 ) ++cnt;
        }

        if( cnt == 0 ) return false;
    }

    // precompute list lengths
    //
    for( size_t i( 0 ); i < qdata_sz_; ++i ) {
        if( qdata_[i].count == 0 ) continue;

        // process high half-word
        //
        {
            TRLI rli( qdata_[i].ext>>SHIFT );
            THalfWord w( 0 );

            int n_words( 0 );

            while( rli.Next( w ) ) {
                if( lut_.serial[w] != serial_ ) {
                    lut_.serial[w] = serial_;
                    lut_.len_l[w] = lut_.len_r[w] = 0;
                    lut_.idx_l[w] = lut_.idx_r[w] = 0xFFFF;
                }

                if( (lut_.idx_l[w]&0xFFF) != (i<<4) ) {
                    ++lut_.len_l[w];
                    lut_.idx_l[w] = 0xF000 + (i<<4);
                }

                ++n_words;
            }
        }

        // process low half-word
        //
        {
            TRLI rli( qdata_[i].ext&MASK );
            THalfWord w( 0 );

            int n_words( 0 );

            while( rli.Next( w ) ) {
                if( lut_.serial[w] != serial_ ) {
                    lut_.serial[w] = serial_;
                    lut_.len_l[w] = lut_.len_r[w] = 0;
                    lut_.idx_l[w] = lut_.idx_r[w] = 0xFFFF;
                }

                if( (lut_.idx_r[w]&0xFFF) != (i<<4) ) {
                    ++lut_.len_r[w];
                    lut_.idx_r[w] = 0xF000 + (i<<4);
                }

                ++n_words;
            }
        }
    }

    // generate offsets and index lists
    //
    common::Uint2 total_len( 0 );

    for( size_t i( 0 ); i < qdata_sz_; ++i ) {
        if( qdata_[i].count == 0 ) continue;

        // process high half-word
        //
        TRLI rli( qdata_[i].ext>>SHIFT );
        THalfWord w( 0 );

        while( rli.Next( w ) ) {
            common::Uint2 & idx( lut_.idx_l[w] );
            common::Uint1 & len( lut_.len_l[w] );

            if( (idx&0xC000) != 0 ) {
                idx = total_len;
                total_len += len;
                len = 0;
            }

            if( len == 0 || q_ext_indices_[idx + len - 1] != i ) {
                q_ext_indices_[idx + (len++)] = i;
            }
        }
    }

    for( size_t i( 0 ); i < qdata_sz_; ++i ) {
        if( qdata_[i].count == 0 ) continue;

        // process low half-word
        //
        TRLI rli( qdata_[i].ext&MASK );
        THalfWord w( 0 );

        while( rli.Next( w ) ) {
            common::Uint2 & idx( lut_.idx_r[w] );
            common::Uint1 & len( lut_.len_r[w] );

            if( (idx&0xC000) != 0 ) {
                idx = total_len;
                total_len += len;
                len = 0;
            }

            if( len == 0 || q_ext_indices_[idx + len - 1] != i ) {
                q_ext_indices_[idx + (len++)] = i;
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
void CBadNMerFilter< t_sdata, t_qdata >::Init( 
        const t_sdata * sdata, const t_qdata * qdata, 
        size_t sdata_sz, size_t qdata_sz, TSeqSize ext_len, int filter_n_err )
{
    SRPRISM_ASSERT( sdata != 0 );
    SRPRISM_ASSERT( qdata != 0 );
    SRPRISM_ASSERT( sdata_sz != 0 );
    SRPRISM_ASSERT( qdata_sz != 0 );
    SRPRISM_ASSERT( qdata_sz <= MAX_QDATA_SIZE );

    sdata_ = sdata;
    qdata_ = qdata;
    sdata_sz_ = sdata_sz;
    qdata_sz_ = qdata_sz;
    ext_len_ = ext_len;
    ext_mask_ = MaskFront< TWord >( ext_len_*LBITS );
    filter_n_err_ = filter_n_err;

    if( ++serial_ == 0 ) { // reinit relevant data
        lut_.serial.assign( SLUTable::LUT_SIZE_1, 0 );
        ++serial_;
    }

    SetState( 
            filter_n_err_ >= (int)ext_len ? 
                (int)SState::ALL : (int)SState::EXACT );

    rnd_map_.clear();
    for( Uint4 i( 0 ); i <= sdata_sz; ++i ) rnd_map_.push_back( i );

    if( randomize_ ) {
        for( Uint4 i( 0 ); i < sdata_sz - 1; ++i ) {
            // Uint4 idx( random()%(sdata_sz - i) );
            Uint4 idx( rand()%(sdata_sz - i) );
            std::swap( rnd_map_[idx], rnd_map_[sdata_sz - i - 1] );
        }
    }
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::NextAll( void )
{
    if( s_ext_idx_ == (ssize_t)sdata_sz_ ) return -1;

    for( size_t i( 0 ); i < qdata_sz_; ++i ) {
        if( qdata_[i].count > 0 ) match_set_.q_ext_set.push_back( i );
    }

    return 0;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::NextExact( void )
{
    if( s_ext_idx_ == (ssize_t)sdata_sz_ ) return -1;
    if( s_ext_idx_ == 0 ) curr_q_ext_ = 0;
    s_extension_ = ((sdata_[match_set_.s_ext].Extension())&ext_mask_);

    while( curr_q_ext_ < qdata_sz_ && qdata_[curr_q_ext_].ext < s_extension_ ) {
        ++curr_q_ext_;
    }

    size_t cqe( curr_q_ext_ );

    while( cqe < qdata_sz_ && 
                qdata_[cqe].ext == s_extension_ ) {
        if( qdata_[cqe].count > 0 ) {
            match_set_.q_ext_set.push_back( cqe );
        }

        ++cqe;
    }

    return 0;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::NextExactRand( void )
{
    if( s_ext_idx_ == (ssize_t)sdata_sz_ ) return -1;
    s_extension_ = ((sdata_[match_set_.s_ext].Extension())&ext_mask_);
    const t_qdata * e( qdata_ + qdata_sz_ ), * r( 0 );
    t_qdata key; key.ext = s_extension_;

    if( (r = std::lower_bound( qdata_, e, key, SQCompare() )) != e ) {
        for( ; r != e && r->ext == s_extension_; ++r ) {
            if( r->count > 0 ) match_set_.q_ext_set.push_back( r - qdata_ );
        }
    }

    return 0;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
template< int n_err > 
inline void CBadNMerFilter< t_sdata, t_qdata >::ProcessLUTEntry( 
        TQExtIdx q_idx )
{
    TWord q_extension( qdata_[q_idx].ext );

    if( dup_checks_[q_idx] != match_set_.s_ext ) {
        dup_checks_[q_idx] = match_set_.s_ext;
        ++stat_;

        if( CFastAlignCheck< TWord, n_err >()( 
                    s_extension_, q_extension, (ext_mask_<<2) ) ) {
            match_set_.q_ext_set.push_back( q_idx );
        }
    }
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
template< int n_err > 
inline void CBadNMerFilter< t_sdata, t_qdata >::ProcessLUTRow( 
        TQExtIdx * start, Uint1 & len )
{
    TQExtIdx * end( start + len );

    while( start < end ) {
        if( qdata_[*start].count == 0 ) {
            --end;
            --len;
            std::swap( *start, *end );
        }
        else {
            ProcessLUTEntry< n_err >( *start );
            ++start;
        }
    }
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
template< int n_err > 
inline void CBadNMerFilter< t_sdata, t_qdata >::NextPriv( void )
{
    static const size_t SHIFT = BYTEBITS*sizeof( THalfWord );
    static const TWord MASK = SBitFieldTraits< TWord, SHIFT >::MASK;

    s_extension_ = ((sdata_[match_set_.s_ext].Extension())&ext_mask_);
    THalfWord lhsw( s_extension_>>SHIFT ), rhsw( s_extension_&MASK ), w( 0 );
        
    if( GetCombinedLenLeft< n_err >( lhsw ) < 
            GetCombinedLenRight< n_err >( rhsw ) ) {
        SRemoveLettersIterator< THalfWord, n_err > RLI( lhsw );

        while( RLI.Next( w ) ) {
            if( lut_.serial[w] == serial_ ) {
                ProcessLUTRow< n_err >( 
                        &q_ext_indices_[0] + lut_.idx_l[w], lut_.len_l[w] );
            }
        }
    }
    else {
        SRemoveLettersIterator< THalfWord, n_err > RLI( rhsw );

        while( RLI.Next( w ) ) {
            if( lut_.serial[w] == serial_ ) {
                ProcessLUTRow< n_err >( 
                        &q_ext_indices_[0] + lut_.idx_r[w], lut_.len_r[w] );
            }
        }
    }
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::NextOneErr( void )
{
    if( s_ext_idx_ == (ssize_t)sdata_sz_ ) return -1;
    NextPriv< 1 >();
    return 1;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::NextTwoErr( void )
{
    if( s_ext_idx_ == (ssize_t)sdata_sz_ ) return -1;
    NextPriv< 2 >();
    return 2;
}

//------------------------------------------------------------------------------
template< typename t_sdata, typename t_qdata >
inline int CBadNMerFilter< t_sdata, t_qdata >::Next( void )
{
    if( state_ == SState::DONE ) return -1;

    if( state_ >= SState::ONE_ERR && s_ext_idx_ == -1 ) {
        if( !SetUpLUT( state_ ) ) { state_ = SState::DONE; return -1; }
    }

    int res( -1 );
    match_set_.q_ext_set.clear();
    match_set_.s_ext = rnd_map_[++s_ext_idx_];

    switch( state_ ) {
        case SState::ALL: 
            res = NextAll(); 
            if( res < 0 ) state_ = SState::DONE;
            break;

        case SState::EXACT:
            res = randomize_ ? NextExactRand() : NextExact();

            if( res < 0 ) {
                if( filter_n_err_ > 0 ) SetState( SState::ONE_ERR );
                else state_ = SState::DONE;
            }

            break;

        case SState::ONE_ERR:
            res = NextOneErr();
            
            if( res < 0 ) { 
                if( filter_n_err_ > 1 ) SetState( SState::TWO_ERR );
                else state_ = SState::DONE;
            }

            break;

        case SState::TWO_ERR: 
            res = NextTwoErr(); 
            if( res < 0 ) state_ = SState::DONE;
            break;

        default: SRPRISM_ASSERT( false );
    }

    return res;
}

END_NS( srprism )
END_STD_SCOPES

#endif

