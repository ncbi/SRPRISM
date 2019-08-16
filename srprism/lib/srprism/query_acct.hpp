/*  $Id: query_acct.hpp 351764 2012-02-01 14:07:34Z morgulis $
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
 * File Description: query rank/level accounting.
 *
 */

#ifndef __SRPRISM_QUERY_ACCT_HPP__
#define __SRPRISM_QUERY_ACCT_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/memmgr.hpp>
#include <srprism/result.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/memmgr.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//##############################################################################
//------------------------------------------------------------------------------
struct CQueryAcct_Base
{
    virtual ~CQueryAcct_Base() {}
};

//##############################################################################
//------------------------------------------------------------------------------
template< typename t_scoring >
class CQueryAcct_TBase : public CQueryAcct_Base
{
    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode DOMAIN_ERR = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == DOMAIN_ERR ) return "domain error";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException );
        };

        typedef t_scoring TScoring;
        typedef typename t_scoring::TQueryData TQueryData;

        static size_t EstimateBytesPerQuery( 
                size_t max_n_res, int n_err, bool paired );

        CQueryAcct_TBase( 
                CMemoryManager & mem_mgr, 
                size_t sz, size_t max_n_res, size_t n_dup, int n_err );
        virtual ~CQueryAcct_TBase();

        void SetResLim( size_t res_lim );

        bool AddResult( 
                TQNum qn, TSeqSize align_len, 
                int n_err, int n_gap, int n_del, int n_gopen );

        bool HasPairedResult( TQNum qn ) const;

        bool AddPairedResult(
                TQNum * qs, TQNum * qe, bool & adjust, TQNum dup_idx,
                TSeqSize alen_1, 
                int n_err_1, int n_gap_1, int n_del_1, int n_gopen_1,
                TSeqSize alen_2, 
                int n_err_2, int n_gap_2, int n_del_2, int n_gopen_2 );

        TQNum AdjustMinRankInfo( TQNum qn, TQNum dup_idx, TQNum curr );

        template< bool paired >
        bool BestLevelFull( TQNum qn ) const
        { return query_info_[qn].template BestLevelFull< paired >( res_lim_ ); }

        template< bool paired > int MaxErr( TQNum qn ) const
        { return query_info_[qn].template MaxErr< paired >(); }

        int GroupMaxErr( TQNum dup_idx ) const
        { return t_scoring::template MaxErr< true >( dup_data_[dup_idx] ); }

        bool HasHigherRank( TQNum qn, const CResult & r ) const
        { return query_info_[qn].HasHigherRank( r ); }

        bool DelResult( TQNum qn, const CResult & r )
        { return query_info_[qn].DelResult( r ); }

        bool AddResult( TQNum qn, const CResult & r );

        void Update( TQNum qn, size_t limit )
        { query_info_[qn].Update( limit ); }

        size_t NRes( TQNum qn ) const 
        { return query_info_[qn].NRes(); }

        void DupScoringData( TQNum dst, TQNum src )
        { query_info_[dst] = query_info_[src]; }

    private:

        typedef typename t_scoring::TRank TRank;

        CMemoryManager & mem_mgr_;
        size_t sz_;
        int n_err_;

        TQueryData * query_info_;
        TRank * dup_data_;
        size_t res_lim_;
};

//------------------------------------------------------------------------------
template< typename t_scoring >
inline size_t CQueryAcct_TBase< t_scoring >::EstimateBytesPerQuery( 
        size_t max_n_res, int n_err, bool paired )
{
    size_t sz( sizeof( TQueryData ) );
    SRPRISM_ASSERT( sz > 0 );
    if( paired ) sz += 1 + sizeof( TQNum )/2;
    sz = ((sz - 1)/4 + 1)*4; // add padding to make sz multiple of 4
    return sz;
}

//------------------------------------------------------------------------------
template< typename t_scoring >
CQueryAcct_TBase< t_scoring >::CQueryAcct_TBase(
        CMemoryManager & mem_mgr, 
        size_t sz, size_t max_n_res, size_t n_dup, int n_err )
    : mem_mgr_( mem_mgr ), sz_( sz ), n_err_( n_err ),
        query_info_( 0 ), dup_data_( 0 ), res_lim_( 0 )
{
    query_info_ = (TQueryData *)mem_mgr_.Allocate( sz_*sizeof( TQueryData ) );
    std::fill( 
            (char *)query_info_, 
            (char *)query_info_ + sz_*sizeof( TQueryData ), 
            0 );

    if( n_dup > 0 ) {
        dup_data_ = (TRank *)mem_mgr_.Allocate( n_dup*sizeof( TRank ) );
        std::fill( 
                (char *)dup_data_, 
                (char *)dup_data_ + n_dup*sizeof( TRank ), 
                0 );
    }
}

//------------------------------------------------------------------------------
template< typename t_scoring >
CQueryAcct_TBase< t_scoring >::~CQueryAcct_TBase()
{
    if( dup_data_ != 0 ) mem_mgr_.Free( (void *)dup_data_ );
    mem_mgr_.Free( (void *)query_info_ );
}

//------------------------------------------------------------------------------
template< typename t_scoring >
void CQueryAcct_TBase< t_scoring >::SetResLim( size_t res_lim )
{
    if( res_lim > t_scoring::MAX_N_RES ) {
        M_THROW( CException, DOMAIN_ERR,
                 "scoring system supports at most " << 
                 t_scoring::MAX_N_RES << " results; " <<
                 res_lim << " requested" );
    }

    res_lim_ = res_lim;
}

//------------------------------------------------------------------------------
template< typename t_scoring >
inline bool CQueryAcct_TBase< t_scoring >::AddResult( 
        TQNum qn, TSeqSize align_len, 
        int n_err, int n_gap, int n_del, int n_gopen )
{
    return query_info_[qn].AddResult( 
            res_lim_, align_len, n_err, n_gap, n_del, n_gopen );
}

//------------------------------------------------------------------------------
template< typename t_scoring >
inline bool CQueryAcct_TBase< t_scoring >::HasPairedResult( TQNum qn ) const
{ return query_info_[qn].HasPairedResult(); }

//------------------------------------------------------------------------------
template< typename t_scoring >
inline bool CQueryAcct_TBase< t_scoring >::AddPairedResult( 
        TQNum * qs, TQNum * qe, bool & adjust, TQNum dup_idx, 
        TSeqSize alen_1, int n_err_1, int n_gap_1, int n_del_1, int n_gopen_1, 
        TSeqSize alen_2, int n_err_2, int n_gap_2, int n_del_2, int n_gopen_2 )
{
    bool has_min_rank( false );
    TQNum qn( *qs );

    if( adjust ) {
        has_min_rank = query_info_[qn].HasEqualRank( dup_data_[dup_idx] );
    }

    if( !query_info_[qn].AddPairedResult( 
                res_lim_, adjust,
                alen_1, n_err_1, n_gap_1, n_del_1, n_gopen_1,
                alen_2, n_err_2, n_gap_2, n_del_2, n_gopen_2 ) ) {
        return false;
    }

    for( TQNum * qi( qs + 1 ); qi != qe; ++qi ) {
        query_info_[*qi] = query_info_[qn];
    }

    adjust = (adjust && has_min_rank);
    return true;
}

//------------------------------------------------------------------------------
template< typename t_scoring >
TQNum CQueryAcct_TBase< t_scoring >::AdjustMinRankInfo( 
        TQNum qn, TQNum dup_idx, TQNum curr )
{
    TQueryData & qd( query_info_[qn] );
    TRank & dd( dup_data_[dup_idx] );
    if( curr == 0 ) { dd = qd.GetRank(); return 1; }

    if( !t_scoring::HasPairedResult( dd ) ) {
        if( !HasPairedResult( qn ) ) return curr + 1;
        else return curr;
    }

    if( qd.HasHigherRank( dd ) ) return curr;
    else if( qd.HasEqualRank( dd ) ) return curr + 1;
    else { dd = qd.GetRank(); return 1; }
}

//------------------------------------------------------------------------------
template< typename t_scoring >
bool CQueryAcct_TBase< t_scoring >::AddResult( TQNum qn, const CResult & r )
{
    if( r.Paired() ) {
        bool dummy;

        return query_info_[qn].AddPairedResult(
                res_lim_, dummy,
                r.GetAlignLen( 0 ), 
                r.NErr( 0 ), r.GetNId( 0 ), r.GetNDel( 0 ), r.GetNGOpen( 0 ),
                r.GetAlignLen( 1 ), 
                r.NErr( 1 ), r.GetNId( 1 ), r.GetNDel( 1 ), r.GetNGOpen( 1 ) );
    }
    else {
        return query_info_[qn].AddResult(
                res_lim_,
                r.GetAlignLen( 0 ), 
                r.NErr( 0 ), r.GetNId( 0 ), r.GetNDel( 0 ), r.GetNGOpen( 0 ) );
    }
}

//##############################################################################
//------------------------------------------------------------------------------
template< typename t_scoring >
class CQueryAcct : public CQueryAcct_TBase< t_scoring >
{
    public:

        CQueryAcct( 
                CMemoryManager & mem_mgr, 
                size_t sz, size_t max_n_res, size_t n_dup, int n_err )
            : CQueryAcct_TBase< t_scoring >( 
                    mem_mgr, sz, max_n_res, n_dup, n_err )
        {
        }
};

END_NS( srprism )
END_STD_SCOPES

#endif

