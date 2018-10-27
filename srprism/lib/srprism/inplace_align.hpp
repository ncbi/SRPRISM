/*  $Id: inplace_align.hpp 573050 2018-10-22 16:56:26Z morgulis $
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
 * File Description: in place alignment for paired searches
 *
 */

#ifndef __PDOP_INPLACE_ALIGN_HPP__
#define __PDOP_INPLACE_ALIGN_HPP__

#include <utility>
#include <vector>
#include <algorithm>

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <common/bits.hpp>
#include <seq/seqdef.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/align.hpp>
#include <srprism/seqstore.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/bits.hpp>
#include <../src/internal/align_toolbox/srprism/lib/seq/seqdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/align.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/stat.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqiter.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CInPlaceAlign
{
    private:

        typedef CSeqStore::TPos TPos;

    public:

        struct Seg
        {
            Seg( TPos p, TSeqSize l, bool lft ) 
                : pos( p ), len( l ), left( lft ) 
            {}

            TPos pos;
            TSeqSize len;
            bool left;
        };

        typedef std::vector< CHit > TResults;
        typedef std::vector< Seg > TSegs;

        CInPlaceAlign( 
                const CSeqStore & ss, const TSegs & segs, 
                TPos anchor_start, TPos anchor_end,
                CExtensionSpaceAllocator & ma, size_t res_lim,
                size_t * n_aligns, size_t * n_ualigns,
                size_t * n_inplace, size_t * n_inplace_aligns )
            : ss_( ss ), segs_( segs ), 
              anchor_start_( anchor_start ), anchor_end_( anchor_end ),
              ma_( ma ), res_lim_( res_lim ), 
              n_aligns_( n_aligns ), n_ualigns_( n_ualigns ),
              n_inplace_( n_inplace ), n_inplace_aligns_( n_inplace_aligns )
        {}

        TPos GetAnchorStart() const { return anchor_start_; }
        TPos GetAnchorEnd() const { return anchor_end_; }

    private:

        bool CheckHit( const CHit & h, const Seg & s ) const
        {
            if( h.Strand() == STRAND_FW ) {
                if( s.left && h.Anchor() > anchor_start_ ) return false;
                if( !s.left && h.Anchor() < anchor_start_ ) return false;
                return (h.Anchor() >= s.pos && 
                        h.Anchor() + h.AnchorSubjLen() <= s.pos + s.len);
            }
            else {
                TPos a( h.Anchor() - h.AnchorSubjLen() );
                if( s.left && a > anchor_start_ ) return false;
                if( !s.left && a < anchor_start_ ) return false;
                return (a >= s.pos && h.Anchor() <= s.pos + s.len);
            }
        }

        template< bool partial_align >
        bool CheckSeed( 
                CHit & hit,
                TSeqSize sa_start, TSeqSize sa_end,
                TSeqSize seed_qoff, TPos seed_soff,
                TSeqSize seed_slen, int n_err )
        {
            ++*n_inplace_aligns_;
            return hit.Extend< partial_align >( 
                    ss_, ma_, sa_start, sa_end, 
                    seed_qoff, seed_soff, seed_slen, 
                    n_err, 0, n_aligns_, n_ualigns_ );
        }

        bool UpdateResults( 
                int & n_err, TResults & results, const CHit & res )
        {
            int res_n_err( res.FullNErr() );

            if( res_n_err < n_err ) {
                results.clear();
                results.push_back( res );
                n_err = res_n_err;
            }
            else {
                bool dup( false );

                for( TResults::reverse_iterator i( results.rbegin() );
                        i != results.rend(); ++i ) {
                    if( i->Strand() == res.Strand() ) {
                        if( res.Strand() == seq::STRAND_FW ) {
                            dup = (i->Anchor() == res.Anchor() ||
                                   i->Anchor() + i->AnchorSubjLen() == 
                                        res.Anchor() + res.AnchorSubjLen());
                        }
                        else {
                            dup = (i->Anchor() == res.Anchor() ||
                                   i->Anchor() - i->AnchorSubjLen() == 
                                        res.Anchor() - res.AnchorSubjLen());
                        }

                        if( dup ) {
                            if( res.FullNId() < i->FullNId() ) *i = res;
                            break;
                        }
                    }
                }

                if( !dup ) results.push_back( res );
            }

            search_n_err_ = n_err;

            if( results.size() == res_lim_ ) {
                if( n_err == 0 ) return true;
                else search_n_err_ = n_err - 1;
            }

            return false;
        }

        template< bool partial_align >
        void SearchLong( 
                const CQueryData & q, TSeqSize sa_start, TSeqSize sa_end,
                int n_err, T_IPAM ipam, TResults & res )
        {
            static const size_t WORD_LETTERS = 
                sizeof( TWord )*seq::SCodingTraits< 
                    SEQDATA_CODING >::PACK_FACTOR;
            static const size_t WORD_LOG_LETTERS = 
                common::SBinLog< WORD_LETTERS >::VALUE;
            static const TSeqSize WORD_LOG_MASK = 
                common::SBitFieldTraits< TSeqSize, WORD_LOG_LETTERS >::MASK;

            int max_seed_err( q.GetSeedNErr() );
            search_n_err_ = n_err;
            int seed_n_err( std::min( n_err, max_seed_err ) );
            int n_hashes( seed_n_err + 1 );
            TWord pat[MAX_N_HASHES];
            bool amb[MAX_N_HASHES];
            bool sa_ambig( q.AmbigSA() );
            TSeqSize sa_off( q.GetSAOffset() );

            {
                CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si( 
                        q.Data(), sa_off, sa_off + n_hashes*HASH_LEN );

                for( int i( 0 ); i < n_hashes; ++i, si += HASH_LEN ) {
                    pat[i] = (*si).first;
                }
            }


            if( sa_ambig ) {
                TSeqSize off( sa_off );

                for( int i( 0 ); i < n_hashes; ++i, off += HASH_LEN ) {
                    amb[i] = (q.NHashAmbigs( off ) > 0);
                }
            }

            for( TSegs::const_iterator i( segs_.begin() );
                    i != segs_.end(); ++i ) {
                if( i->len + search_n_err_ < q.Len() ) continue;
                T_IPAM lcl_ipam( 0 );

                if( i->left ) {
                    if( (ipam&IPAM_LEFT_ENABLED) == 0 ) continue;
                    else lcl_ipam = (ipam&IPAM_LEFT_ENABLED);
                }
                else {
                    if( (ipam&IPAM_RIGHT_ENABLED) == 0 ) continue;
                    else lcl_ipam = (ipam&IPAM_RIGHT_ENABLED);
                }

                typedef std::pair< const TWord *, TSeqSize > TSSPtr;
                TPos p( i->pos );
                TSSPtr fw_ptr( ss_.FwDataPtr( p ) ),
                       rv_ptr( ss_.RvDataPtr( p + WORD_LETTERS ) );

                for( TSeqSize s( 0 ); 
                        s < i->len - WORD_LETTERS; ++s ) {
                    TWord fw_subj_word( 
                            seq::GetWord< SEQDATA_CODING, TWord >(
                                fw_ptr.first, fw_ptr.second ) ),
                          rv_subj_word(
                            seq::GetWord< SEQDATA_CODING, TWord >(
                                rv_ptr.first, rv_ptr.second ) );

                    // check forward seeds
                    //
                    if( (lcl_ipam&IPAM_FW_ENABLED) != 0 ) {
                        TSeqSize off( sa_off );

                        for( int j( 0 ); j < n_hashes; ++j, off += HASH_LEN ) {
                            if( sa_ambig && amb[j] ) continue;

                            if( fw_subj_word == pat[j] ) {
                                CHit r( q, seq::STRAND_FW );

                                if( CheckSeed< partial_align >( 
                                            r, sa_start, sa_end, off, p + s, 
                                            HASH_LEN, search_n_err_ ) ) {
                                    if( CheckHit( r, *i ) ) {
                                        if( UpdateResults( n_err, res, r ) ) {
                                            return;
                                        }

                                        seed_n_err = std::min(
                                                max_seed_err, search_n_err_ );
                                        n_hashes = seed_n_err + 1;
                                    }
                                }
                            }
                        }
                    }

                    // check reverse seeds
                    //
                    if( (lcl_ipam&IPAM_RV_ENABLED) != 0 ) {
                        TSeqSize off( sa_off );

                        for( int j( 0 ); j < n_hashes; ++j, off += HASH_LEN ) {
                            if( sa_ambig && amb[j] ) continue;

                            if( rv_subj_word == pat[j] ) {
                                CHit r( q, seq::STRAND_RV );

                                if( CheckSeed< partial_align >( 
                                            r, sa_start, sa_end, off, p + s, 
                                            HASH_LEN, search_n_err_ ) ) {
                                    if( CheckHit( r, *i ) ) {
                                        if( UpdateResults( n_err, res, r ) ) {
                                            return;
                                        }

                                        seed_n_err = std::min(
                                                max_seed_err, search_n_err_ );
                                        n_hashes = seed_n_err + 1;
                                    }
                                }
                            }
                        }
                    }

                    fw_ptr.second = ((fw_ptr.second + 1)&WORD_LOG_MASK);
                    rv_ptr.second = ((rv_ptr.second - 1)&WORD_LOG_MASK);
                    fw_ptr.first += 
                        ((WORD_LETTERS - fw_ptr.second)>>WORD_LOG_LETTERS);
                    rv_ptr.first -= ((rv_ptr.second + 1)>>WORD_LOG_LETTERS);
                }
            }
        }

    public:

        template< bool partial_align >
        void Search( 
                const CQueryData & q, TSeqSize sa_start, TSeqSize sa_end,
                int n_err, T_IPAM ipam, TResults & res )
        {
            ++*n_inplace_;

            if( q.IsSALong() && q.GetNHashes() > 2 ) {
                SearchLong< partial_align >( 
                        q, sa_start, sa_end, n_err, ipam, res );
                return;
            }

            static const size_t WORD_LETTERS = 
                sizeof( TWord )*seq::SCodingTraits< 
                    SEQDATA_CODING >::PACK_FACTOR;
            static const size_t HALF_WORD_LETTERS = (WORD_LETTERS>>1);

            static const size_t HALF_WORD_BITS = 
                (common::SIntTraits< TWord >::BITS>>1);
            static const size_t WORD_LOG_LETTERS = 
                common::SBinLog< WORD_LETTERS >::VALUE;
            static const TSeqSize WORD_LOG_MASK = 
                common::SBitFieldTraits< TSeqSize, WORD_LOG_LETTERS >::MASK;
            static const TWord HALF_WORD_MASK = 
                SBitFieldTraits< TWord, HALF_WORD_BITS >::MASK;

            static const TWord MASKS[3] = 
                { 0xFFFFFFFFULL, 0xFFFFFFFFULL, 0xFFFFULL };
            static const size_t SHIFTS[3] = { 0, 0, HALF_WORD_BITS };
            static const TSeqSize LENGTHS[3] = { 
                WORD_LETTERS, WORD_LETTERS, HALF_WORD_LETTERS };

            static const size_t N_MATCHES_0 = 1;
            static const size_t N_MATCHES_1 = 2;
            static const size_t N_MATCHES_2 = 4;
            static const size_t N_MATCHES[3] = { 
                N_MATCHES_0, N_MATCHES_1, N_MATCHES_2 };

            static const TSeqSize OFFSETS_0[N_MATCHES_0] = { 0 };
            static const TSeqSize OFFSETS_1[N_MATCHES_1] = { 0, WORD_LETTERS };
            static const TSeqSize OFFSETS_2[N_MATCHES_2] = { 
                0, HALF_WORD_LETTERS, 
                2*HALF_WORD_LETTERS, 3*HALF_WORD_LETTERS };
            static const TSeqSize * OFFSETS[3] = { 
                OFFSETS_0, OFFSETS_1, OFFSETS_2 };

            size_t short_adj( q.Len() < MIN_MED_QUERY_LEN ? 1 : 0 );

            int max_seed_err( q.GetSeedNErr() );
            search_n_err_ = n_err;
            int seed_n_err( std::min( n_err, max_seed_err ) );

            // prepare query prefix patterns
            //
            TWord pat_0[N_MATCHES_0],
                  pat_1[N_MATCHES_1],
                  pat_2[N_MATCHES_2];
            TWord * pat_p[3] = { pat_0, pat_1, pat_2 };
            TSeqSize sa_off( q.GetSAOffset() );
            CSeqFwIterator< SEQDATA_CODING, TWord, TWord > si(
                    q.Data(), sa_off, sa_off + HASH_LEN + HASH_LEN );

            pat_0[0] = (*si).first;
            pat_1[0] = pat_0[0];
            si += HASH_LEN;
            pat_1[1] = short_adj > 0 ? 0 : (*si).first;
            pat_2[0] = (pat_1[0]>>HALF_WORD_BITS);
            pat_2[1] = (pat_1[0]&HALF_WORD_MASK);
            pat_2[2] = (pat_1[1]>>HALF_WORD_BITS);
            pat_2[3] = (pat_1[1]&HALF_WORD_MASK);

            bool amb_0[N_MATCHES_0] = {false},
                 amb_1[N_MATCHES_1] = {false, false},
                 amb_2[N_MATCHES_2] = {false, false, false};
            bool * amb_p[3] = { amb_0, amb_1, amb_2 };
            bool sa_ambig( q.AmbigSA() );

            if( sa_ambig ) {
                amb_2[0] = (q.NAmbigs( sa_off, sa_off + HASH_LEN/2 ) > 0);
                amb_2[1] = (q.NAmbigs( 
                            sa_off + HASH_LEN/2, sa_off + HASH_LEN) > 0);
                amb_0[0] = amb_1[0] = (amb_2[0] || amb_2[1]);

                if( !short_adj ) {
                    amb_2[2] = (q.NAmbigs( 
                                sa_off + HASH_LEN, 
                                sa_off + HASH_LEN + HASH_LEN/2 ) > 0);
                    amb_2[3] = (q.NAmbigs(
                                sa_off + HASH_LEN + HASH_LEN/2,
                                sa_off + HASH_LEN + HASH_LEN ) > 0);
                    amb_1[1] = (amb_2[2] || amb_2[3]);
                }
            }

            CQueryData q_copy( q );
            q_copy.SetType1( SErrType::E );
            q_copy.SetType2( SErrType::E );

            // search seeds
            //
            for( TSegs::const_iterator i( segs_.begin() );
                    i != segs_.end(); ++i ) {
                if( i->len + search_n_err_ < q.Len() ) continue;
                T_IPAM lcl_ipam( 0 );

                if( i->left ) {
                    if( (ipam&IPAM_LEFT_ENABLED) == 0 ) continue;
                    else lcl_ipam = (ipam&IPAM_LEFT_ENABLED);
                }
                else {
                    if( (ipam&IPAM_RIGHT_ENABLED) == 0 ) continue;
                    else lcl_ipam = (ipam&IPAM_RIGHT_ENABLED);
                }

                typedef std::pair< const TWord *, TSeqSize > TSSPtr;
                TPos p( i->pos );
                TSSPtr fw_ptr( ss_.FwDataPtr( p ) ),
                       rv_ptr( ss_.RvDataPtr( p + WORD_LETTERS ) );
                TWord mask( MASKS[seed_n_err + short_adj] );
                size_t shift( SHIFTS[seed_n_err + short_adj] );
                TSeqSize length( LENGTHS[seed_n_err + short_adj] );
                TWord * pat( pat_p[seed_n_err + short_adj] );
                bool * amb( amb_p[seed_n_err + short_adj] );
                const TSeqSize * offsets( OFFSETS[seed_n_err + short_adj] );
                size_t n_matches( N_MATCHES[seed_n_err] );

                for( TSeqSize s( 0 ); 
                        s < i->len - HALF_WORD_LETTERS; ++s ) {
                    TWord fw_subj_word( 
                            seq::GetWord< SEQDATA_CODING, TWord >(
                                fw_ptr.first, fw_ptr.second ) ),
                          rv_subj_word(
                            seq::GetWord< SEQDATA_CODING, TWord >(
                                rv_ptr.first, rv_ptr.second ) );

                    // check forward seeds
                    //
                    if( (lcl_ipam&IPAM_FW_ENABLED) != 0 ) {
                        for( size_t j( 0 ); j < n_matches; ++j ) {
                            if( sa_ambig && amb[j] ) continue;

                            if( (fw_subj_word>>shift) == pat[j] ) {
                                CHit r( q_copy, seq::STRAND_FW );

                                if( CheckSeed< partial_align >( 
                                            r, sa_start, sa_end,
                                            sa_off + offsets[j], p + s, 
                                            length, search_n_err_ ) ) {
                                    if( CheckHit( r, *i ) ) {
                                        if( UpdateResults( n_err, res, r ) ) {
                                            return;
                                        }

                                        seed_n_err = std::min( 
                                                max_seed_err, search_n_err_ );
                                        mask = MASKS[seed_n_err + short_adj];
                                        shift = SHIFTS[seed_n_err + short_adj];
                                        length = 
                                            LENGTHS[seed_n_err + short_adj];
                                        pat = pat_p[seed_n_err + short_adj];
                                        amb = amb_p[seed_n_err + short_adj];
                                        offsets = 
                                            OFFSETS[seed_n_err + short_adj];
                                        n_matches = N_MATCHES[seed_n_err];
                                    }
                                }
                            }
                        }
                    }

                    // check reverse seeds
                    //
                    if( (lcl_ipam&IPAM_RV_ENABLED) != 0 ) {
                        for( size_t j( 0 ); j < n_matches; ++j ) {
                            if( sa_ambig && amb[j] ) continue;

                            if( (rv_subj_word&mask) == pat[j] ) {
                                CHit r( q_copy, seq::STRAND_RV );

                                if( CheckSeed< partial_align >( 
                                            r, sa_start, sa_end,
                                            sa_off + offsets[j], p + s, 
                                            length, search_n_err_ ) ) {
                                    if( CheckHit( r, *i ) ) {
                                        if( UpdateResults( n_err, res, r ) ) {
                                            return;
                                        }

                                        seed_n_err = std::min( 
                                                max_seed_err, search_n_err_ );
                                        mask = MASKS[seed_n_err + short_adj];
                                        shift = SHIFTS[seed_n_err + short_adj];
                                        length = 
                                            LENGTHS[seed_n_err + short_adj];
                                        pat = pat_p[seed_n_err + short_adj];
                                        amb = amb_p[seed_n_err + short_adj];
                                        offsets = 
                                            OFFSETS[seed_n_err + short_adj];
                                        n_matches = N_MATCHES[seed_n_err];
                                    }
                                }
                            }
                        }
                    }

                    fw_ptr.second = ((fw_ptr.second + 1)&WORD_LOG_MASK);
                    rv_ptr.second = ((rv_ptr.second - 1)&WORD_LOG_MASK);
                    fw_ptr.first += 
                        ((WORD_LETTERS - fw_ptr.second)>>WORD_LOG_LETTERS);
                    rv_ptr.first -= ((rv_ptr.second + 1)>>WORD_LOG_LETTERS);
                }
            }
        }

    private:

        const CSeqStore & ss_;
        const TSegs & segs_;
        TPos anchor_start_;
        TPos anchor_end_;
        CExtensionSpaceAllocator & ma_;
        size_t res_lim_;
        int search_n_err_;
        size_t * n_aligns_;
        size_t * n_ualigns_;
        size_t * n_inplace_;
        size_t * n_inplace_aligns_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

