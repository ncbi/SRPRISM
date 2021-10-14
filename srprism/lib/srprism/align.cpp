/*  $Id: align.cpp 639115 2021-10-13 15:24:22Z morgulis $
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

#include <ncbi_pch.hpp>

#include <iterator>
#include <algorithm>

#include "align.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
inline TWord CHit::CutWord( TDWord d, const TSeqSize & o )
{
    TSeqSize ob( WBITS - (o<<LSHIFT) );
    return (TWord)((d>>ob)&LWORD_MASK);
}

//------------------------------------------------------------------------------
inline void CHit::ExtendLeftExact( 
        const TWord * q, const TWord * qm, const TWord * s, const TWord * sm,
        TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
        TDWord dq, TDWord dm, TDWord ds, TDWord dsm, TSeqSize & n_matched )
{
    if( s_lim > 0 ) {
        TWord qw, mw, sw, smw;
        TSeqSize l( WLETTERS );

        do {
            qw = CutWord( dq, qwoff );
            mw = CutWord( dm, qwoff );
            sw = CutWord( ds, swoff );
            smw = CutWord( dsm, swoff );
            qw ^= sw; qw |= mw; qw |= smw;
            l = (common::FirstSetBit_Right( qw )>>LSHIFT);
            n_matched += l;
            s_lim -= l;
            if( s_lim <= 0 || l != WLETTERS ) break;
            --q; --qm; --s; --sm;
            dq >>= WBITS; dm >>= WBITS; ds >>= WBITS; dsm >>= WBITS;
            dq += (((TDWord)(*q))<<WBITS);
            dm += (((TDWord)(*qm))<<WBITS);
            ds += (((TDWord)(*s))<<WBITS);
            dsm += (((TDWord)(*sm))<<WBITS);
        }
        while( true );

        if( s_lim < 0 ) n_matched += s_lim;
    }
}

//------------------------------------------------------------------------------
inline void CHit::ExtendRightExact( 
        const TWord * q, const TWord * qm, const TWord * s, const TWord * sm,
        TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
        TDWord dq, TDWord dm, TDWord ds, TDWord dsm, TSeqSize & n_matched )
{
    if( s_lim > 0 ) {
        TWord qw, mw, sw, smw;
        TSeqSize l( WLETTERS );
        ++q; ++qm; ++s; ++sm;

        do {
            qw = CutWord( dq, qwoff );
            mw = CutWord( dm, qwoff );
            sw = CutWord( ds, swoff );
            smw = CutWord( dsm, swoff );
            qw ^= sw; qw |= mw; qw |= smw;
            l = (common::FirstSetBit_Left( qw )>>LSHIFT);
            n_matched += l;
            s_lim -= l;
            if( s_lim <= 0 || l != WLETTERS ) break;
            ++q; ++qm; ++s; ++sm;
            dq <<= WBITS; dm <<= WBITS; ds <<= WBITS; dsm <<= WBITS;
            dq += *q; dm += *qm; ds += *s; dsm += *sm;
        }
        while( true );

        if( s_lim < 0 ) n_matched += s_lim;
    }
}

//------------------------------------------------------------------------------
#define VISITED(idx) (mat[idx].serial == serial)

#define VISIT(idx) SRPRISM_ASSERT( idx < mat_sz ); \
    if( !VISITED( idx ) ) { \
        mat[idx].penalty = 0xFFFF; \
        mat[idx].range = 0; \
        mat[idx].serial = serial; \
    }

#define ADJUST_NODE SRPRISM_ASSERT( idx < mat_sz ); \
                    SRPRISM_ASSERT( curr_idx < mat_sz ); \
                    SRPRISM_ASSERT( last_shoot < mat_sz ); \
                    if( !VISITED( idx ) || mat[idx].penalty > penalty ) { \
                        if( !VISITED( idx ) ) { \
                            shoots[last_shoot++] = idx; \
                            mat[idx].range = 0; \
                            mat[idx].serial = serial; \
                        } \
                        mat[idx].penalty = penalty; \
                        mat[curr_idx].range = range; \
                    }

//------------------------------------------------------------------------------

template<> struct CHit::SExtendDirTraits< -1 >
{
    static void ExtendExact(
        const TWord * q, const TWord * qm, const TWord * s, const TWord * sm,
        TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
        TDWord dq, TDWord dm, TDWord ds, TDWord dsm, TSeqSize & n_matched )
    {
        ExtendLeftExact( q, qm, s, sm, qwoff, swoff, s_lim,
                         dq, dm, ds, dsm, n_matched );
    }

    static void AdjustDataPtr( 
            TSeqSize off_ini, Uint2 off_curr, 
            Uint2 & res_adj, TSeqSize & res_off )
    {
        off_ini -= off_curr;
        res_adj = ((off_ini)>>WSHIFT);
        res_off = (off_ini&WMASK);
    }
};

template<> struct CHit::SExtendDirTraits< 1 >
{
    static void ExtendExact(
        const TWord * q, const TWord * qm, const TWord * s, const TWord * sm,
        TSeqSize qwoff, TSeqSize swoff, common::Sint4 s_lim, 
        TDWord dq, TDWord dm, TDWord ds, TDWord dsm, TSeqSize & n_matched )
    {
        ExtendRightExact( q, qm, s, sm, qwoff, swoff, s_lim,
                          dq, dm, ds, dsm, n_matched );
    }

    static void AdjustDataPtr( 
            TSeqSize off_ini, Uint2 off_curr, 
            Uint2 & res_adj, TSeqSize & res_off )
    {
        off_ini += off_curr;
        res_adj = (off_ini>>WSHIFT);
        res_off = (off_ini&WMASK);
    }
};

//------------------------------------------------------------------------------
inline void CHit::STraceBack::TrimAndUpdate( 
        Uint4 curr_idx, TSeqSize q_lim, TSeqSize q_min, TSeqSize q_max, 
        TSeqSize range, int n_err )
{
    const SMatrixEntry * c( mat + curr_idx ), 
                       * e;
    int c_n_err( c->penalty>>8 );
    TSeqSize qpos;

    if( c_n_err > n_err ) {
        qpos = curr_idx/col_len;
        if( qpos + 1 <= q_min ) return;
    }
    else qpos = q_lim - 1;

    int score;

    while( true ) {
        if( 1 + qpos < q_min ) break;
        score = 1 + qpos - c_n_err;

        if( score >= best_score ) {
            best_idx = c - mat - (range<<col_shift);
            best_score = score;
            flank = q_max - best_idx/col_len - range - 1;
        }

        if( (int)range > c_n_err ||
                (c - (range<<col_shift))->penalty == 0 ) {
            break;
        }

        qpos -= range;
        c -= (range<<col_shift);
        e = c - col_len;

        if( e > mat && e->serial == serial ) {
            e -= (e->range<<col_shift);

            if( e->serial == serial && e->penalty + ERR_STEP == c->penalty ) {
                c -= col_len;
                --qpos;
                range = c->range;
                --c_n_err;
                continue;
            }
        }

        e = c + 1;

        if( e->serial == serial ) {
            e -= (e->range<<col_shift);

            if( e->serial == serial && 
                    e->penalty + ERR_STEP + 1 == c->penalty ) {
                ++c;
                range = c->range;
                --c_n_err;
                continue;
            }
        }

        e = c - col_len - 1;

        if( e > mat && e->serial == serial ) {
            e -= (e->range<<col_shift);

            if( e->serial == serial && 
                    e->penalty + ERR_STEP + 1 == c->penalty ) {
                --qpos;
                c -= col_len + 1;
                range = c->range;
                --c_n_err;
                continue;
            }
        }

        SRPRISM_ASSERT( false );
    }
}

//------------------------------------------------------------------------------
template< int dir, bool partial_align >
inline void CHit::STraceBack::ExtendNonExact(
        const TWord * q, const TWord * m, const TWord * s, const TWord * sm,
        TSeqSize qwoff, TSeqSize swoff, 
        TSeqSize q_min, TSeqSize q_max, TSeqSize q_lim, TSeqSize s_lim, 
        common::Uint1 n_err, common::Uint1 n_id, 
        CExtensionSpaceAllocator & ma )
{
    SRPRISM_ASSERT( q_lim > 1 );
    SRPRISM_ASSERT( q_min <= q_lim );
    flank = 0;

    // matrix dimensions
    //
    col_shift = BinLog( 2*(n_id + 1) ) + 1;
    col_len = (1<<col_shift);
    TSeqSize col_mask( col_len - 1 );

// #ifndef NDEBUG
#if defined( _DEBUG ) || !defined( NDEBUG )
    SRPRISM_ASSERT( col_len <= 255 );
    size_t mat_sz( q_max*col_len );
#endif

    // work space allocation
    //
    Uint4 * shoots;
    serial = ma.GetSerial();

    {
        CExtensionSpaceAllocator::TExtensionSpaceHandle h( 
                ma.Alloc( q_max, col_len ) );
        mat = h.first;
        shoots = h.second;
    }

    // other initialization
    //
    best_score = 0;
    best_idx = 0;
    start_idx = n_id + 1;
    Uint4 curr_idx( start_idx ), idx;
    Uint4 next_shoot( 0 ), last_shoot( 0 );
    VISIT( curr_idx );
    mat[start_idx].penalty = 0;
    idx = curr_idx + col_len;

    if( s_lim > 1 ) {
        VISIT( idx );
        mat[idx].penalty = ERR_STEP;
    }

    // rest of initialization
    //
    if( s_lim > 1 ) shoots[last_shoot++] = idx;

    ++idx;
    VISIT( idx );
    mat[idx].penalty = ERR_STEP + 1;
    shoots[last_shoot++] = idx;

    if( s_lim > 1 ) {
        idx -= col_len + 2;
        VISIT( idx );
        mat[idx].penalty = ERR_STEP + 1;
        shoots[last_shoot++] = idx;
    }

    Uint2 qadj, sadj, qadj_p, sadj_p, i, j;
    TSeqSize range;
    TSeqSize qo, so;
    SExtendDirTraits< dir >::AdjustDataPtr( qwoff, 0, qadj_p, qo );
    SExtendDirTraits< dir >::AdjustDataPtr( swoff, 0, sadj_p, so );
    TDWord dq( *(q + qadj_p) ), dm( *(m + qadj_p) ),
           ds( *(s + sadj_p) ), dsm( *(sm + sadj_p) );
    dq = (dq<<WBITS) + *(q + qadj_p + 1);
    dm = (dm<<WBITS) + *(m + qadj_p + 1);
    ds = (ds<<WBITS) + *(s + sadj_p + 1);
    dsm = (dsm<<WBITS) + *(sm + sadj_p + 1);

    // main loop
    //
    while( last_shoot > next_shoot ) {
        curr_idx = shoots[next_shoot++];
        i = (curr_idx>>col_shift);
        Uint2 p( (mat[curr_idx].penalty>>8) );

        if( partial_align ) {
            if( i + 1U >= q_lim || p > n_err ) {
                TrimAndUpdate( curr_idx, q_lim, q_min, q_max, 0, n_err );
                continue;
            }
        }
        else {
            if( i + 1U >= q_lim || p > n_err ) {
                best_idx = curr_idx;
                break;
            }
        }

        j = i + n_id + 1 - (curr_idx&col_mask);
        range = 0;

        SExtendDirTraits< dir >::AdjustDataPtr( qwoff, i + 1, qadj, qo );
        SExtendDirTraits< dir >::AdjustDataPtr( swoff, j + 1, sadj, so );

        if( qadj != qadj_p ) {
            dq = *(q + qadj); dq = (dq<<WBITS) + *(q + qadj + 1);
            dm = *(m + qadj); dm = (dm<<WBITS) + *(m + qadj + 1);
            qadj_p = qadj;
        }

        if( sadj != sadj_p ) {
            ds  = *(s + sadj);  ds  = (ds<<WBITS)  + *(s + sadj + 1);
            dsm = *(sm + sadj); dsm = (dsm<<WBITS) + *(sm + sadj + 1);
            sadj_p = sadj;
        }

        SExtendDirTraits< dir >::ExtendExact( 
                q + qadj, m + qadj, s + sadj, sm + sadj, qo, so, 
                (Sint4)(s_lim - j - 1), dq, dm, ds, dsm, range );

        if( i + range + 1 >= q_lim ) {
            if( partial_align ) {
                Uint2 penalty( mat[curr_idx].penalty );
                curr_idx += (range<<col_shift);
                mat[curr_idx].penalty = penalty;
                TrimAndUpdate( curr_idx, q_lim, q_min, q_max, range, n_err );
            }
            else best_idx = curr_idx;

            break;
        }

        Uint2 penalty( mat[curr_idx].penalty + ERR_STEP );
        curr_idx += (range<<col_shift);
        VISIT( curr_idx );
        idx = curr_idx + col_len;
        i += range + 1; j += range + 1;
        if( j < s_lim ) { ADJUST_NODE; }
        ++penalty;

        if( i + 1 <= n_id + 1 + j ) {
            ++idx; 
            ADJUST_NODE; 
            --idx;
        }

        if( j + 1 <= n_id + 1 + i && j < s_lim ) {
            idx -= col_len + 1;
            ADJUST_NODE;
        }
    }
}

//------------------------------------------------------------------------------
#undef ADJUST_NODE
#undef VISIT
#undef VISITED

//------------------------------------------------------------------------------
template< bool partial_align >
bool CHit::Extend( 
        const CSeqStore & ss, CExtensionSpaceAllocator & ma, 
        TSeqSize sa_start, TSeqSize sa_end,
        TSeqSize seed_qoff, TPos seed_soff, TSeqSize seed_slen, 
        int n_err, CBreakSegs * bsegs, size_t * stat, size_t * ustat )
{
    ++*stat;
    int seed_n_err( std::min( n_err, qdata_.GetSeedNErr() ) - qdata_.NErr() );

    // initializations and declarations
    //
    TSeqSize seed_qlen( seed_slen + qdata_.NIns() - qdata_.NDel() );
    TSeqSize s_fw_lim( ss.FwTailLen( seed_soff ) - seed_slen ),
             s_rv_lim( ss.RvTailLen( seed_soff ) ),
             q_fw_lim( qdata_.Len() - seed_qoff - seed_qlen ),
             q_rv_lim( seed_qoff );
    const TWord * q( qdata_.Data() ), * m( qdata_.Mask() ), * s, * sm;
    TSeqSize qo_l( seed_qoff ), qo_r( seed_qoff + seed_qlen ),
             so_l, so_r;
    int dir;

    // strand dependent initialization
    //
    if( strand_ == STRAND_FW ) {
        std::pair< const TWord *, TSeqSize > sd( ss.FwDataPtr( seed_soff ) );
        s = sd.first; so_l = sd.second; so_r = sd.second + seed_slen;
        sm = ss.FwMaskPtr( seed_soff );
        dir = -1;
    }
    else {
        std::swap( s_fw_lim, s_rv_lim );
        seed_soff += seed_slen;
        std::pair< const TWord *, TSeqSize > sd( ss.RvDataPtr( seed_soff ) );
        s = sd.first; so_l = sd.second; so_r = sd.second + seed_slen;
        sm = ss.RvMaskPtr( seed_soff );
        dir = 1;
    }

    // adjusting subject start pointer and offset
    //
    {
        TSeqSize l( std::min( s_rv_lim, q_rv_lim + n_err ) );
        l = (l>>WSHIFT) + 1;
        s -= l;
        sm -= l;
        l <<= WSHIFT;
        so_l += l; so_r += l;
    }

    align_len_ = left_flank_ = 0;

    // initial left exact extension
    //
    n_err -= qdata_.NErr();
    bool check_left( false );

    if( qo_l > 0 ) {
        Uint2 qa( (qo_l)>>WSHIFT ), sa( (so_l)>>WSHIFT );
        const TWord * qq( q + qa ), * mm( m + qa ),
                    * ss( s + sa ), * ssm( sm + sa );
        TDWord dq( *(qq - 1) ), dm( *(mm - 1) ),
               ds( *(ss - 1) ), dsm( *(ssm - 1) );
        dq = (dq<<WBITS) + *qq;
        dm = (dm<<WBITS) + *mm;
        ds = (ds<<WBITS) + *ss;
        dsm = (dsm<<WBITS) + *ssm;
        TSeqSize exact_extension_left( 0 );
        ExtendLeftExact(
                qq - 1, mm - 1, ss - 1, ssm - 1,
                (qo_l&WMASK), (so_l&WMASK), s_rv_lim,
                dq, dm, ds, dsm, exact_extension_left );

        if( !partial_align && exact_extension_left < seed_qoff && n_err == 0 ) {
            return false;
        }

        check_left = 
            (bsegs != 0 && !bsegs->Empty() && bsegs->Left( 0 ) < seed_qoff );
        
        if( check_left && exact_extension_left > 0 ) {
            bsegs->SetBounds( 
                    seed_qoff - exact_extension_left, seed_qoff + seed_qlen );
            if( !bsegs->Check() ) return false;
        }

        qo_l -= exact_extension_left;
        so_l -= exact_extension_left;
        q_rv_lim -= exact_extension_left;
        s_rv_lim -= exact_extension_left;
    }

    // initial right exact extension
    //
    bool check_right( false );

    if( q_fw_lim > 0 ) {
        Uint2 qa( qo_r>>WSHIFT ), sa( so_r>>WSHIFT );
        const TWord * qq( q + qa ), * mm( m + qa ),
                    * ss( s + sa ), * ssm( sm + sa );
        TDWord dq( *qq ), dm( *mm ), ds( *ss ), dsm( *ssm );
        dq = (dq<<WBITS) + *(qq + 1);
        dm = (dm<<WBITS) + *(mm + 1);
        ds = (ds<<WBITS) + *(ss + 1);
        dsm = (dsm<<WBITS) + *(ssm + 1);
        TSeqSize exact_extension_right( 0 );
        ExtendRightExact(
                qq, mm, ss, ssm, (qo_r&WMASK), (so_r&WMASK), s_fw_lim,
                dq, dm, ds, dsm, exact_extension_right );

        if( !partial_align && exact_extension_right < q_fw_lim && n_err == 0 ) {
            return false;
        }

        check_right = 
            (bsegs != 0 && !bsegs->Empty() && 
             bsegs->Right( bsegs->Size() - 1 ) >= seed_qoff );

        if( check_right && q_fw_lim > 0 ) {
            bsegs->SetBounds( seed_qoff, qo_r + exact_extension_right );
            if( !bsegs->Check() ) return false;
        }

        qo_r += exact_extension_right;
        so_r += exact_extension_right;
        q_fw_lim -= exact_extension_right;
        s_fw_lim -= exact_extension_right;
    }

    left_flank_ = q_rv_lim;
    align_len_ = qdata_.Len() - left_flank_ - q_fw_lim;
    err_data_.reserve( n_err + qdata_.NErr() );
    anchor_ = seed_soff;
    sa_end = std::min( sa_end, qdata_.Len() );

    bool left_exact( sa_start == 0 ),
         right_exact( sa_end == qdata_.Len() );
    int left_sa_n_err( 0 );
    TSeqSize left_q_lim( 
            q_rv_lim + 1 > sa_start ? q_rv_lim + 1 - sa_start : 0 );

    // initial left extension to the boundary of the seeding area
    //
    if( q_rv_lim > sa_start ) {
        ma.Clean();
        STraceBack tb;
        ++*ustat;
        tb.ExtendNonExact< -1, false >(
                q - 1, m - 1, s - 1, sm - 1, qo_l + 1, so_l + 1,
                0, q_rv_lim + 1, left_q_lim, s_rv_lim + 1,
                seed_n_err, seed_n_err, ma );
        if( (left_sa_n_err = tb.NErr()) > seed_n_err ) return false;
        if( left_exact ) { n_err -= tb.NErr(); n_id_ += tb.NId(); }
        if( check_left ) bsegs->SetBounds( 0, seed_qoff + seed_qlen );

        if( check_left || left_exact ) {
            SErrRec e;
            tb.Reset();

            while( tb.NextError( e ) ) {
                e.qpos = q_rv_lim - e.qpos;
                bool d( e.err_type == SErrType::D );

                if( d ) {
                    --e.qpos;
                    if( left_exact ) { anchor_ += dir; ++n_d_; }
                }
                else if( left_exact && e.err_type == SErrType::I ) {
                    anchor_ -= dir;
                }

                if( check_left ) bsegs->SetCurr( e.qpos, d, true );
                if( left_exact ) err_data_.push_back( e );
            }

            if( check_left && !bsegs->Check() ) return false;
        }
    }

    TSeqSize right_q_lim( 
            q_fw_lim + 1 > qdata_.Len() - sa_end ?  
                q_fw_lim - (qdata_.Len() - sa_end) + 1 : 0 );
    int right_sa_n_err( 0 );
    TErrors right_err_data;

    // initial right extension to the boundary of the seeding area
    //
    if( q_fw_lim > qdata_.Len() - sa_end ) {
        ma.Clean();
        STraceBack tb;
        int n_err_lcl( seed_n_err - left_sa_n_err );
        ++*ustat;
        tb.ExtendNonExact< 1, false >(
                q, m, s, sm, qo_r - 1, so_r - 1,
                0, q_fw_lim + 1, right_q_lim, s_fw_lim + 1,
                n_err_lcl, n_err_lcl, ma );
        if( (right_sa_n_err = tb.NErr()) > n_err_lcl ) return false;
        if( check_right ) bsegs->SetBounds( seed_qoff, qdata_.Len() );

        if( right_exact ) {
            n_id_ += tb.NId();
            right_err_data.resize( right_sa_n_err );
        }

        if( right_exact || check_right ) {
            SErrRec e;
            tb.Reset();
            size_t err_idx( right_err_data.size() );

            while( tb.NextError( e ) ) {
                e.qpos += qo_r - 1;
                bool d( e.err_type == SErrType::D );

                if( check_right ) {
                    bsegs->SetCurr( e.qpos + (d ? 1 : 0), d, false );
                }

                if( right_exact ) {
                    right_err_data[--err_idx] = e;
                    if( d ) ++n_d_;
                }
            }

            if( check_right && !bsegs->Check() ) return false;
        }
    }

    n_err -= right_sa_n_err;

    // left extension to the end of the query
    //
    if( q_rv_lim > 0 && !left_exact ) {
        ma.Clean();
        STraceBack tb;
        ++*ustat;
        tb.ExtendNonExact< -1, partial_align >(
                q - 1, m - 1, s - 1, sm - 1, qo_l + 1, so_l + 1,
                left_q_lim, q_rv_lim + 1, q_rv_lim + 1, s_rv_lim + 1,
                n_err, n_err, ma );

        if( partial_align ) {
            SRPRISM_ASSERT( tb.NErr() <= n_err );
            if( tb.best_idx == 0 ) return false;
        }
        else if( tb.NErr() > n_err ) return false;

        n_err -= tb.NErr();
        n_id_ += tb.NId();
        SErrRec e;
        tb.Reset();

        while( tb.NextError( e ) ) {
            e.qpos = q_rv_lim - e.qpos;

            if( e.err_type == SErrType::D ) {
                anchor_ += dir;
                --e.qpos;
                ++n_d_;
            }
            else if( e.err_type == SErrType::I ) anchor_ -= dir;

            err_data_.push_back( e );
        }

        left_flank_ = tb.flank;
    }
    else left_flank_ = 0;

    // append internal seed error data
    //
    if( qdata_.NErr() > 0 ) {
        SErrRec e = { (Uint2)(seed_qoff + qdata_.GetPos1()), qdata_.Type1() };
        if( e.err_type == SErrType::D ) --e.qpos;
        if( !qdata_.RightBUDir() ) e.qpos += qdata_.NIns();
        err_data_.push_back( e );
    }

    if( qdata_.NErr() > 1 ) {
        SErrRec e = { (Uint2)(seed_qoff + qdata_.GetPos2()), qdata_.Type2() };
        if( !qdata_.RightBUDir() ) e.qpos += qdata_.NIns();
        err_data_.push_back( e );
    }

    if( right_exact ) {
        std::copy( 
                right_err_data.begin(), right_err_data.end(),
                std::back_inserter( err_data_ ) );
    }

    n_err += right_sa_n_err;

    // right extension to the end of the query
    //
    if( q_fw_lim > 0 && !right_exact ) {
        if( !partial_align && n_err == 0 ) return false;
        ma.Clean();
        STraceBack tb;
        ++*ustat;
        tb.ExtendNonExact< 1, partial_align >( 
                q, m, s, sm, qo_r - 1, so_r - 1,
                right_q_lim, q_fw_lim + 1, q_fw_lim + 1, s_fw_lim + 1,
                n_err, n_err, ma );
        
        if( partial_align ) {
            SRPRISM_ASSERT( tb.NErr() <= n_err );
            if( tb.best_idx == 0 ) return false;
        }
        else if( tb.NErr() > n_err ) return false;

        n_id_ += tb.NId();
        err_data_.resize( err_data_.size() + tb.NErr() );
        size_t err_idx( err_data_.size() );
        tb.Reset();

        while( err_idx > 0 && tb.NextError( err_data_[--err_idx] ) ) {
            SErrRec & e( err_data_[err_idx] );
            e.qpos += qo_r - 1;
            bool d( e.err_type == SErrType::D );
            if( d ) ++n_d_;
        }

        align_len_ = qdata_.Len() - left_flank_ - tb.flank + n_d_;
    }
    else align_len_ = qdata_.Len() - left_flank_ + n_d_;

    anchor_ += dir*(seed_qoff - left_flank_);

    if( !partial_align ) {
        left_flank_ = 0;
        align_len_ = qdata_.Len() + n_d_;
    }

    return true;
}

// explicitly instantiate CHit::Extend()
//
template bool CHit::Extend< false >( 
        const CSeqStore &, CExtensionSpaceAllocator &, 
        TSeqSize, TSeqSize,
        TSeqSize, TPos, TSeqSize, int, CBreakSegs *, size_t *, size_t * );

template bool CHit::Extend< true >( 
        const CSeqStore &, CExtensionSpaceAllocator &, 
        TSeqSize, TSeqSize,
        TSeqSize, TPos, TSeqSize, int, CBreakSegs *, size_t *, size_t * );

END_NS( srprism )
END_STD_SCOPES

