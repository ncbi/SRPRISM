/*  $Id: scoring.hpp 590234 2019-07-25 16:28:25Z morgulis $
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
 * File Description: scoring system interface and implementations
 *
 */

#ifndef __SRPRISM_SCORING_HPP__
#define __SRPRISM_SCORING_HPP__

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <common/exception.hpp>
#include <srprism/srprismdef.hpp>
#include <srprism/result.hpp>
#include <srprism/align.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/exception.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/result.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/align.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
/*
   REQUIREMENTS FOR SCORING SYSTEM CLASSES

   Ranks and Levels: 

    - only results with the highest rank are kept;
    - results of the lower level are preferred.

   typedef ... TRank;

        Integer type identifying a rank of a result.

   static const size_t MAX_N_ERR;

        Maximum number of errors supported by the scoring system.

   static const size_t MAX_N_RES;

        Maximum number of results supported by the scoring system.

   TQueryData - holds fixed size per-query information

   static bool HaveEqualRanks( const CResult & l, const CResult & r );

        Returns true iff two results have the same rank.

   static bool HaveEqualLevels( const CResult & l, const CResult & r );

        Returns true iff two results have equal levels (assumes
        that they have equal ranks).

   static bool CompareByRank( const CResult & l, const CResult & r );

        Returns true iff l is higher ranked than r.

   static bool CompareByLevel( const CResult & l, const CResult & r );

        Returns true iff l has lower level than r
        (assumes that l and r have equal ranks).

   static bool CompareByLevelAtCol( 
        const CResult & l, const CResult & r, size_t col_idx );

        Returns true iff l has lower level than r at column col_idx
        (assumes that l and r have equal ranks).

   template< bool paired > static int MaxErr( TRank rank );

        Return maximum number of errors that a query with a given
        rank can have.

   bool HasPairedResult( TRank r );

        True iff a result with given rank corresponds to a paired alignment.

   ========================================
   Interface to TQueryData

        bool AddResult( 
            size_t limit, TSeqSize align_len, 
            int n_err, int n_gap, int n_del, int n_gopen );

        bool AddPairedResult(
            size_t limit, bool & adjust,
            TSeqSize alen_1, int n_err_1, int n_gap_1, int n_del_1, int n_gopen_1,
            TSeqSize alen_2, int n_err_2, int n_gap_2, int n_del_2, int n_gopen_2 );

            Register single or paired alignment with the query.
            Return true on success; false if request fails due to limit
            being reached. 'adjust' should be set true if the
            rank of the query inreased due to addition; false otherwise

        bool DelResult( const CResult & r );
            
            De-register a result from the query data.

        TRank GetRank( void ) const;

            Return an integer representation of query rank. This
            should only be used by acct classes.

        bool HasPairedResult( void ) const;

            Returns true if there is a paired result registered for
            the query.

        bool HasEqualRank( const TQueryData & rhs ) const;
        bool HasEqualRank( TRank rhs ) const;

            Returns true if this query has the same rank as rhs.

        bool HasHigherRank( const TQueryData & rhs ) const;
        bool HasHigherRank( const CResult & rhs ) const;
        bool HasHigherRank( TRank rhs ) const;

            Returns true if this query has higher rank than rhs.

        TRank GetRank( void ) const;

            Get rank value for the query.

        void Update( size_t limit );

            Adjust query information based on the new limit on the
            number of results.

        size_t NRes( void ) const;

            Return the number of results registered for this query.

        template< bool paired > int MaxErr( void ) const;

            Maximum number of errors a query can have in a result.

        template< bool paired > bool BestLevelFull( size_t lim ) const;

            Returns true if the lim resuts are found for the current
            best rank.
*/

//##############################################################################
//
// Exceptions of the scoring systems.
//
class ScoringSystemException : public common::CException
{
    private:

        typedef common::CException TBase;

    public:

        static const TErrorCode DOMAIN_ERROR = 1;

        virtual const std::string ErrorMessage( TErrorCode code ) const
        {
            if( code == DOMAIN_ERROR ) return "domain error";
            else return TBase::ErrorMessage( code );
        }

        M_EXCEPT_CTOR( ScoringSystemException );
};

//##############################################################################
//------------------------------------------------------------------------------
/*
    SCORING SYSTEM FOR DEFAULT GLOBAL SEARCH

    Single alignment rank is 127 - #errors
    Paired-end alignemt rank is 255 - max( e1, e2 ), where e1, e2 are number
        of errors in the first and the second mate alignments respectively.

    Single alignment level = # of indels
    Paired-end alignment level is min( l1, l2 ), where l1 and l2 are number
        of indels in the first and the second mate alignments respectively.
*/
class CDefaultScoringSystem
{
    public:

        static const size_t MAX_N_ERR = 16;
        typedef common::Uint2 TRank;

    private:

        typedef common::Uint1 TLevelCount;

        struct SSUUtil
        {
            typedef common::Uint1 TSubUnit;
            typedef common::SIntTraits< TSubUnit > TTraits;

            static const size_t BITS = TTraits::BITS;
            static const size_t MAX = TTraits::MAX;
            static const size_t MASK = MAX;
        };

        template< common::Uint1 base_rank >
        struct SMakeRank
        {
            static const TRank VALUE = 
                (base_rank<<SSUUtil::BITS) + SSUUtil::MAX;
        };

        static const TRank MAX_RANK = ((SSUUtil::MAX)<<SSUUtil::BITS);
        static const TRank MAX_UNPAIRED_RANK = SMakeRank< 127 >::VALUE;
        static const TRank BASE_RANK_UNIT = (1ULL<<SSUUtil::BITS);

        static TRank NErr2Rank( int n_err )
        { return MAX_UNPAIRED_RANK - n_err*BASE_RANK_UNIT; }

        static TRank NErr2Rank( int n_err_1, int n_err_2 )
        {
            return (MAX_RANK - std::max( n_err_1, n_err_2 )*BASE_RANK_UNIT) +
                   (BASE_RANK_UNIT - 1) - (n_err_1 + n_err_2);
        }

    private:

        static TRank NErr2BestRank( int n_err ) 
        { return NErr2Rank( n_err, 0 ); }

    public:

        struct CQueryData
        {
            bool AddResult( 
                    size_t limit, TSeqSize align_len, 
                    int n_err, int n_gap, int n_del, int n_gopen );

            bool HasPairedResult( void ) const;

            bool HasEqualRank( const CQueryData & r ) const
            { return rank_ == r.rank_; }

            bool HasEqualRank( TRank rank ) const { return rank_ == rank; }

            TRank GetRank( void ) const { return rank_; }

            bool AddPairedResult(
                    size_t limit, bool & adjust,
                    TSeqSize, int n_err_1, int n_gap_1, int n_del_1, int,
                    TSeqSize, int n_err_2, int n_gap_2, int n_del_2, int );

            bool HasHigherRank( const CQueryData & r ) const
            { return rank_ > r.rank_; }

            bool HasHigherRank( const CResult & r ) const;

            bool HasHigherRank( TRank r ) const { return rank_ > r; }

            template< bool paired > bool BestLevelFull( size_t lim ) const;

            template< bool paired > int MaxErr( void ) const;

            bool DelResult( const CResult & r );

            void Update( size_t limit );

            size_t NRes( void ) const;

            private:

                void ClearCounts( void )
                { std::fill( counts_, counts_ + MAX_N_ERR, 0 ); }

                bool IncrementCounts( int n_gap, size_t limit );

                bool AddResultPriv( 
                        size_t limit, bool & adjust, TRank rank, int level );

                TRank rank_;
                common::Uint1 n_res_;
                TLevelCount counts_[MAX_N_ERR];
        };

        static TRank ResultRank( const CResult & r );
        static int ResultLevel( const CResult & r );

    public:

        static const size_t MAX_N_RES = common::SIntTraits< TLevelCount >::MAX;

        typedef CQueryData TQueryData;

        static bool HaveEqualRanks( const CResult & l, const CResult & r );
        static bool HaveEqualLevels( const CResult & l, const CResult & r );
        static bool CompareByRank( const CResult & l, const CResult & r );
        static bool CompareByLevel( const CResult & l, const CResult & r );
        static bool CompareByLevelAtCol( 
                const CResult & l, const CResult & r, size_t col_idx );

        template< bool paired > static int MaxErr( TRank rank );
        static bool HasPairedResult( TRank r );
};

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::IncrementCounts( 
        int n_gap, size_t limit )
{
    if( n_res_ < limit ) { ++n_res_; ++counts_[n_gap]; }
    else {
        int i( n_gap );
        for( ; i > 0 && counts_[i-1] == 0; --i );

        if( i > n_gap + 1 ) {
            --counts_[i-1];
            ++counts_[n_gap];
        }
        else return false;
    }

    SRPRISM_ASSERT( n_res_ <= limit );
    return true;
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::AddResultPriv( 
        size_t limit, bool & adjust, TRank rank, int level )
{
    adjust = false;

    if( rank > rank_ ) {
        rank_ = rank;
        ClearCounts();
        counts_[level] = 1;
        n_res_ = 1;
        adjust = true;
        return true;
    }
    else if( rank == rank_ ) return IncrementCounts( level, limit );
    else return false;
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::AddResult( 
        size_t limit, TSeqSize , int n_err, int n_gap, int n_del, int )
{
    TRank rank( NErr2Rank( n_err ) );
    bool dummy;
    return AddResultPriv( limit, dummy, rank, n_gap );
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::AddPairedResult( 
        size_t limit, bool & adjust,
        TSeqSize, int n_err_1, int n_gap_1, int n_del_1, int, 
        TSeqSize, int n_err_2, int n_gap_2, int n_del_2, int )
{
    TRank rank( NErr2Rank( n_err_1, n_err_2 ) );
    int level( std::min( n_gap_1, n_gap_2 ) );
    return AddResultPriv( limit, adjust, rank, level );
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::HasPairedResult( void ) const
{ return (rank_ > MAX_UNPAIRED_RANK); }

//------------------------------------------------------------------------------
template<> inline int 
CDefaultScoringSystem::CQueryData::MaxErr< false >( void ) const
{ return (MAX_UNPAIRED_RANK - (rank_&~SSUUtil::MASK))>>SSUUtil::BITS; }

//------------------------------------------------------------------------------
template<> inline int 
CDefaultScoringSystem::CQueryData::MaxErr< true >( void ) const
{ 
    return rank_ > MAX_UNPAIRED_RANK ? 
                        ((MAX_RANK - (rank_&~SSUUtil::MASK))>>SSUUtil::BITS) : 
                        (MAX_UNPAIRED_RANK>>SSUUtil::BITS);
}

//------------------------------------------------------------------------------
template<> inline bool
CDefaultScoringSystem::CQueryData::BestLevelFull< false >( size_t lim ) const
{ return counts_[0] == lim; }

//------------------------------------------------------------------------------
template<> inline bool
CDefaultScoringSystem::CQueryData::BestLevelFull< true >( size_t lim ) const
{
    return (rank_ == NErr2BestRank( MaxErr< true >() )) && counts_[0] == lim; 
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::HasHigherRank( 
        const CResult & r ) const
{
    TRank rrank( r.Paired() ? NErr2Rank( r.NErr( 0 ), r.NErr( 1 ) )
                            : NErr2Rank( r.NErr( 0 ) ) );
    return rank_ > rrank;
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CQueryData::DelResult( const CResult & r )
{
    int level( r.Paired() ? std::min( r.GetNId( 0 ), r.GetNId( 1 ) )
                          : r.GetNId( 0 ) );
    if( counts_[level] == 0 ) return false;
    --counts_[level];
    --n_res_;
    return true;
}

//------------------------------------------------------------------------------
inline void CDefaultScoringSystem::CQueryData::Update( size_t limit )
{
    size_t s( 0 );

    for( size_t i( 0 ); i < MAX_N_ERR; ++i ) {
        if( s + counts_[i] > limit ) {
            counts_[i] = limit - s;
            s = limit;
        }
        else s += counts_[i];
    }
}

//------------------------------------------------------------------------------
inline size_t CDefaultScoringSystem::CQueryData::NRes( void ) const
{ return n_res_; }

//------------------------------------------------------------------------------
inline CDefaultScoringSystem::TRank 
CDefaultScoringSystem::ResultRank( const CResult & r )
{
    return r.Paired() ? NErr2Rank( r.NErr( 0 ), r.NErr( 1 ) ) 
                      : NErr2Rank( r.NErr( 0 ) );
}

//------------------------------------------------------------------------------
inline int CDefaultScoringSystem::ResultLevel( const CResult & r )
{
    if( r.Paired() ) return std::min( r.GetNId( 0 ), r.GetNId( 1 ) );
    else return r.GetNId( 0 );
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::HaveEqualRanks( 
        const CResult & l, const CResult & r )
{ return ResultRank( l ) == ResultRank( r ); }

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::HaveEqualLevels( 
        const CResult & l, const CResult & r )
{ 
    SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
    return ResultLevel( l ) == ResultLevel( r );
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CompareByRank( 
        const CResult & l, const CResult & r )
{ return ResultRank( l ) > ResultRank( r ); }

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CompareByLevel( 
        const CResult & l, const CResult & r )
{
    SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
    return ResultLevel( l ) < ResultLevel( r );
}

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::CompareByLevelAtCol( 
        const CResult & l, const CResult & r, size_t col_idx )
{
    SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
    return l.GetNId( col_idx ) < r.GetNId( col_idx );
}

//------------------------------------------------------------------------------
template<> inline int CDefaultScoringSystem::MaxErr< true >( TRank rank )
{ 
    return rank > MAX_UNPAIRED_RANK ? 
                    ((MAX_RANK - (rank&~SSUUtil::MASK))>>SSUUtil::BITS) : 
                    (MAX_UNPAIRED_RANK>>SSUUtil::BITS); 
}

template<> inline int CDefaultScoringSystem::MaxErr< false >( TRank rank )
{ return (MAX_UNPAIRED_RANK - (rank&~SSUUtil::MASK))>>SSUUtil::BITS; }

//------------------------------------------------------------------------------
inline bool CDefaultScoringSystem::HasPairedResult( TRank r )
{ return r > MAX_UNPAIRED_RANK; }

//##############################################################################
//------------------------------------------------------------------------------
/*
    SCORING SYSTEM FOR PARTIAL_ALIGNMENT SEARCH

    Rank of an alignment is the number of query positions participating
    in the alignment minus number of errors.

    Rank of paired-end alignment is the sum of ranks of the two mate alignments.

    Single alignment level = # of indels
    Paired-end alignment level is min( l1, l2 ), where l1 and l2 are number
        of indels in the first and the second mate alignments respectively.
*/
class CLenMinusErrScoringSystem
{
    private:

        struct SRank
        {
            SRank( common::Sint2 a  = 0, bool hp = false )
                : rank( a ), has_paired( hp )
            {}

            common::Sint2 rank;
            bool has_paired;

            friend bool operator==( const SRank & l, const SRank & r )
            { return l.rank == r.rank && l.has_paired == r.has_paired; }

            friend bool operator!=( const SRank & l, const SRank & r )
            { return !( l == r ); }

            friend bool operator<( const SRank & l, const SRank & r )
            {
                if( l.has_paired == r.has_paired )
                {
                    return l.rank < r.rank;
                }

                return r.has_paired;
            }

            friend bool operator>( const SRank & l, const SRank & r )
            { return r < l; }
        };

    public:

        static const size_t MAX_N_ERR = 16;
        typedef SRank TRank;

    private:

        typedef common::Uint1 TLevelCount;

        struct CQueryData
        {
            bool AddResult( 
                    size_t limit, TSeqSize align_len, 
                    int n_err, int n_gap, int n_del, int n_gopen );

            bool HasPairedResult( void ) const { return rank_.has_paired; }

            bool HasEqualRank( const CQueryData & r ) const
            { return rank_ == r.rank_; }

            bool HasEqualRank( TRank rank ) const { return rank_ == rank; }

            TRank GetRank( void ) const { return rank_; }

            bool AddPairedResult(
                    size_t limit, bool & adjust,
                    TSeqSize, int n_err_1, int n_gap_1, int n_del_1, int,
                    TSeqSize, int n_err_2, int n_gap_2, int n_del_2, int );

            bool HasHigherRank( const CQueryData & r ) const
            { return rank_ > r.rank_; }

            bool HasHigherRank( const CResult & r ) const
            { return rank_ > ResultRank( r ); }

            bool HasHigherRank( TRank r ) const { return rank_ > r; }

            template< bool paired > bool BestLevelFull( size_t lim ) const
            { return false; }

            template< bool paired > int MaxErr( void ) const { return MAX_N_ERR; }

            bool DelResult( const CResult & r );

            void Update( size_t limit );

            size_t NRes( void ) const { return n_res_; }

            private:

                void ClearCounts( void )
                { std::fill( counts_, counts_ + MAX_N_ERR, 0 ); }

                bool IncrementCounts( int n_gap, size_t limit );

                bool AddResultPriv( 
                        size_t limit, bool & adjust, TRank rank, int level );

                TRank rank_;
                common::Uint1 n_res_;
                TLevelCount counts_[MAX_N_ERR];
        };

        static TRank ResultRank( const CResult & r );

        static int ResultLevel( const CResult & r )
        { 
            return r.Paired() ? std::min( r.GetNId( 0 ), r.GetNId( 1 ) ) 
                              : r.GetNId( 0 );
        }

    public:

        static const size_t MAX_N_RES = common::SIntTraits< TLevelCount >::MAX;

        typedef CQueryData TQueryData;

        static bool HaveEqualRanks( const CResult & l, const CResult & r )
        { return ResultRank( l ) == ResultRank( r ); }

        static bool CompareByRank( const CResult & l, const CResult & r )
        { return ResultRank( l ) > ResultRank( r ); }

        static bool HaveEqualLevels( const CResult & l, const CResult & r )
        { return ResultLevel( l ) == ResultLevel( r ); }

        static bool CompareByLevel( const CResult & l, const CResult & r )
        { return ResultLevel( l ) < ResultLevel( r ); }

        static bool CompareByLevelAtCol( 
                const CResult & l, const CResult & r, size_t col_idx )
        { return l.GetNId( col_idx ) < r.GetNId( col_idx ); }

        template< bool paired > static int MaxErr( TRank rank ) { return MAX_N_ERR; }
        static bool HasPairedResult( TRank r ) { return r.has_paired; }
};

//------------------------------------------------------------------------------
inline bool CLenMinusErrScoringSystem::CQueryData::IncrementCounts( 
        int n_gap, size_t limit )
{
    if( n_res_ < limit ) { ++n_res_; ++counts_[n_gap]; }
    else {
        int i( n_gap );
        for( ; i > 0 && counts_[i-1] == 0; --i );

        if( i > n_gap + 1 ) {
            --counts_[i-1];
            ++counts_[n_gap];
        }
        else return false;
    }

    SRPRISM_ASSERT( n_res_ <= limit );
    return true;
}

//------------------------------------------------------------------------------
inline bool CLenMinusErrScoringSystem::CQueryData::AddResultPriv( 
        size_t limit, bool & adjust, TRank rank, int level )
{
    adjust = false;

    if( rank > rank_ ) {
        rank_ = rank;
        ClearCounts();
        counts_[level] = 1;
        n_res_ = 1;
        adjust = true;
        return true;
    }
    else if( rank == rank_ ) return IncrementCounts( level, limit );
    else return false;
}

//------------------------------------------------------------------------------
inline bool CLenMinusErrScoringSystem::CQueryData::AddResult( 
        size_t limit, TSeqSize align_len, int n_err, int n_gap, int n_del, int )
{
    TRank rank( align_len - n_del - n_err );
    bool dummy;
    return AddResultPriv( limit, dummy, rank, n_gap );
}

//------------------------------------------------------------------------------
inline bool CLenMinusErrScoringSystem::CQueryData::AddPairedResult( 
        size_t limit, bool & adjust,
        TSeqSize align_len_1, int n_err_1, int n_gap_1, int n_del_1, int, 
        TSeqSize align_len_2, int n_err_2, int n_gap_2, int n_del_2, int )
{
    common::Sint2 s1( align_len_1 - n_del_1 - n_err_1 ),
                  s2( align_len_2 - n_del_2 - n_err_2 );
    TRank rank( s1 + s2, true );
    int level( std::min( n_gap_1, n_gap_2 ) );
    return AddResultPriv( limit, adjust, rank, level );
}

//------------------------------------------------------------------------------
inline bool CLenMinusErrScoringSystem::CQueryData::DelResult( 
        const CResult & r )
{
    int level( r.Paired() ? std::min( r.GetNId( 0 ), r.GetNId( 1 ) )
                          : r.GetNId( 0 ) );
    if( counts_[level] == 0 ) return false;
    --counts_[level];
    --n_res_;
    return true;
}

//------------------------------------------------------------------------------
inline void CLenMinusErrScoringSystem::CQueryData::Update( size_t limit )
{
    size_t s( 0 );

    for( size_t i( 0 ); i < MAX_N_ERR; ++i ) {
        if( s + counts_[i] > limit ) {
            counts_[i] = limit - s;
            s = limit;
        }
        else s += counts_[i];
    }
}

//------------------------------------------------------------------------------
inline CLenMinusErrScoringSystem::TRank 
CLenMinusErrScoringSystem::ResultRank( const CResult & r )
{
    TRank res;
    res.rank = r.GetAlignLen( 0 ) - r.GetNDel( 0 ) - r.NErr( 0 );

    if( r.Paired() ) {
        common::Sint2 s( r.GetAlignLen( 1 ) - r.GetNDel( 1 ) - r.NErr( 1 ) );
        res.rank += s;
        res.has_paired = true;
    }

    return res;
}

//##############################################################################
//------------------------------------------------------------------------------
/*
    SCORING SYSTEM FOR GLOBAL SEARCH WITH SUM OF #(MATE ERRORS) FOR PAIRED
    RESULTS

    Single alignment rank is 127 - #errors
    Paired-end alignemt rank is 255 - (e1 + e2), where e1, e2 are number
        of errors in the first and the second mate alignments respectively.

    Single alignment level = # of indels
    Paired-end alignment level is min( l1, l2 ), where l1 and l2 are number
        of indels in the first and the second mate alignments respectively.
*/
class CSumErrScoringSystem
{
    public:

        static const size_t MAX_N_ERR = 16;
        typedef common::Uint1 TRank;
        typedef common::Uint1 TLevelCount;
        static const size_t MAX_N_RES = common::SIntTraits< TLevelCount >::MAX;

    private:

        static const TRank MAX_UNPAIRED_RANK = 127;
        static const TRank MAX_RANK = 255;

        static TRank NErr2Rank( int n_err ) 
        { return MAX_UNPAIRED_RANK - n_err; }

        static TRank NErr2Rank( int n_err_1, int n_err_2 )
        { return MAX_RANK - (n_err_1 + n_err_2 ); }

        static TRank ResultRank( const CResult & r )
        {
            return r.Paired() ? NErr2Rank( r.NErr( 0 ), r.NErr( 1 ) )
                              : NErr2Rank( r.NErr( 0 ) );
        }

        static int ResultLevel( const CResult & r )
        {
            return r.Paired() ? std::min( r.GetNId( 0 ), r.GetNId( 1 ) )
                              : r.GetNId( 0 );
        }

    public:

        struct CQueryData
        {
            bool AddResult( 
                    size_t limit, TSeqSize align_len, 
                    int n_err, int n_gap, int, int n_gopen )
            {
                TRank rank( NErr2Rank( n_err ) );
                bool dummy;
                return AddResultPriv( limit, dummy, rank, n_gap );
            }

            bool AddPairedResult(
                    size_t limit, bool & adjust,
                    TSeqSize, int n_err_1, int n_gap_1, int, int,
                    TSeqSize, int n_err_2, int n_gap_2, int, int )
            {
                TRank rank( NErr2Rank( n_err_1, n_err_2 ) );
                int level( std::min( n_gap_1, n_gap_2 ) );
                return AddResultPriv( limit, adjust, rank, level );
            }

            bool DelResult( const CResult & r )
            {
                int level( ResultLevel( r ) );
                if( counts_[level] == 0 ) return false;
                --counts_[level];
                --n_res_;
                return true;
            }

            bool HasPairedResult( void ) const
            { return (rank_ > MAX_UNPAIRED_RANK); }

            TRank GetRank( void ) const { return rank_; }

            bool HasEqualRank( const CQueryData & r ) const
            { return rank_ == r.rank_; }

            bool HasEqualRank( TRank rank ) const { return rank_ == rank; }

            bool HasHigherRank( const CQueryData & r ) const
            { return rank_ > r.rank_; }

            bool HasHigherRank( TRank r ) const { return rank_ > r; }

            bool HasHigherRank( const CResult & r ) const
            {
                return HasHigherRank(
                        r.Paired() ? NErr2Rank( r.NErr( 0 ), r.NErr( 1 ) )
                                   : NErr2Rank( r.NErr( 0 ) ) );
            }

            void Update( size_t limit )
            {
                size_t s( 0 );

                for( size_t i( 0 ); i < MAX_N_ERR; ++i ) {
                    if( s + counts_[i] > limit ) {
                        counts_[i] = limit - s;
                        s = limit;
                    }
                    else s += counts_[i];
                }
            }

            size_t NRes( void ) const { return n_res_; }

            template< bool paired > int MaxErr( void ) const
            { return CSumErrScoringSystem::MaxErr< paired >( rank_ ); }

            template< bool paired > bool BestLevelFull( size_t lim ) const
            { return counts_[0] == lim; }

            private:

                void ClearCounts( void )
                { std::fill( counts_, counts_ + MAX_N_ERR, 0 ); }

                bool IncrementCounts( int n_gap, size_t limit )
                {
                    if( n_res_ < limit ) { ++n_res_; ++counts_[n_gap]; }
                    else {
                        int i( n_gap );
                        for( ; i > 0 && counts_[i-1] == 0; --i );
                
                        if( i > n_gap + 1 ) {
                            --counts_[i-1];
                            ++counts_[n_gap];
                        }
                        else return false;
                    }

                    SRPRISM_ASSERT( n_res_ <= limit );
                    return true;
                }

                bool AddResultPriv( 
                        size_t limit, bool & adjust, TRank rank, int level )
                {
                    adjust = false;

                    if( rank > rank_ ) {
                        rank_ = rank;
                        ClearCounts();
                        counts_[level] = 1;
                        n_res_ = 1;
                        adjust = true;
                        return true;
                    }
                    else if( rank == rank_ ) {
                        return IncrementCounts( level, limit );
                    }
                    else return false;
                }

                TRank rank_;
                common::Uint1 n_res_;
                common::Uint1 counts_[MAX_N_ERR];
        };

        typedef CQueryData TQueryData;

        static bool HaveEqualRanks( const CResult & l, const CResult & r )
        { return ResultRank( l ) == ResultRank( r ); }

        static bool HaveEqualLevels( const CResult & l, const CResult & r )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return ResultLevel( l ) == ResultLevel( r );
        }

        static bool CompareByRank( const CResult & l, const CResult & r )
        { return ResultRank( l ) > ResultRank( r ); }

        static bool CompareByLevel( const CResult & l, const CResult & r )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return ResultLevel( l ) < ResultLevel( r );
        }

        static bool CompareByLevelAtCol( 
                const CResult & l, const CResult & r, size_t col_idx )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return l.GetNId( col_idx ) < r.GetNId( col_idx );
        }

        template< bool paired > static int MaxErr( TRank rank );

        static bool HasPairedResult( TRank r )
        { return r > MAX_UNPAIRED_RANK; }
};

//------------------------------------------------------------------------------
template<> inline int CSumErrScoringSystem::MaxErr< false >( TRank rank )
{ return MAX_UNPAIRED_RANK - rank; }

template<> inline int CSumErrScoringSystem::MaxErr< true >( TRank rank )
{ return (rank > MAX_UNPAIRED_RANK) ? MAX_RANK - rank : MAX_UNPAIRED_RANK; }

//##############################################################################
/*
   SCORING SYSTEM FOR REPORTING RESULTS WITH AT MOST A GIVEN NUMBER OF ERRORS

   All results have the same rank; level values correspond to the number
   of errors (e.g. mismatches now are not preferred to indels).

   For paired-end alignments the level is the max of the levels of individual
   mate alignments.
*/
//------------------------------------------------------------------------------
class CWeakScoringSystem
{
    public:

        static const size_t MAX_N_ERR = 16;
        typedef common::Uint1 TRank;

    private:

        static const TRank RANK_SINGLE = 0;
        static const TRank RANK_PAIRED = 1;

        typedef common::Uint1 TLevelCount;

    public:

        struct CQueryData
        {
            bool AddResult( 
                    size_t limit, TSeqSize align_len, 
                    int n_err, int , int , int )
            {
                bool dummy;
                return AddResultPriv( limit, dummy, RANK_SINGLE, n_err );
            }

            bool AddPairedResult(
                    size_t limit, bool & adjust,
                    TSeqSize, int n_err_1, int , int , int,
                    TSeqSize, int n_err_2, int , int , int )
            {
                return AddResultPriv( 
                        limit, adjust, RANK_PAIRED,
                        std::min( n_err_1, n_err_2 ) );
            }

            bool HasPairedResult( void ) const 
            { return (rank_ == RANK_PAIRED); }

            bool HasEqualRank( const CQueryData & r ) const
            { return (rank_ == r.rank_); }

            bool HasEqualRank( TRank rank ) const { return (rank_ == rank); }

            TRank GetRank( void ) const { return rank_; }

            bool HasHigherRank( const CQueryData & r ) const
            { return (rank_ > r.rank_); }

            bool HasHigherRank( const CResult & r ) const
            { return rank_ > ResultRank( r ); }

            bool HasHigherRank( TRank r ) const { return (rank_ > r); }

            template< bool paired > bool BestLevelFull( size_t lim ) const;

            template< bool paired > int MaxErr( void ) const
            { return MAX_N_ERR; }

            bool DelResult( const CResult & r )
            {
                int level( ResultLevel( r ) );

                if( counts_[level] == 0 ) return false;
                --counts_[level];
                --n_res_;
                return true;
            }

            void Update( size_t limit )
            {
                size_t s( 0 );

                for( size_t i( 0 ); i < MAX_N_ERR; ++i ) {
                    if( s + counts_[i] > limit ) {
                        counts_[i] = limit - s;
                        s = limit;
                    }
                    else s += counts_[i];
                }
            }

            size_t NRes( void ) const { return n_res_; }

            private:

                void ClearCounts( void )
                { std::fill( counts_, counts_ + MAX_N_ERR, 0 ); }

                bool IncrementCounts( int n_err, size_t limit )
                {
                    if( n_res_ < limit ) { ++n_res_; ++counts_[n_err]; }
                    else {
                        int i( n_err );
                        for( ; i > 0 && counts_[i-1] == 0; --i );

                        if( i > n_err + 1 ) {
                            --counts_[i-1];
                            ++counts_[n_err];
                        }
                        else return false;
                    }

                    SRPRISM_ASSERT( n_res_ <= limit );
                    return true;
                }

                bool AddResultPriv( 
                        size_t limit, bool & adjust, TRank rank, int level )
                {
                    adjust = false;

                    if( rank > rank_ ) {
                        rank_ = rank;
                        ClearCounts();
                        counts_[level] = 1;
                        n_res_ = 1;
                        adjust = true;
                        return true;
                    }
                    else if( rank == rank_ ) {
                        return IncrementCounts( level, limit );
                    }
                    else return false;
                }

                TRank rank_;
                common::Uint1 n_res_;
                TLevelCount counts_[MAX_N_ERR];
        };

        static TRank ResultRank( const CResult & r )
        { return r.Paired() ? (TRank)RANK_PAIRED : (TRank)RANK_SINGLE; }

        static int ResultLevel( const CResult & r )
        { 
            return r.Paired() ? std::min( r.NErr( 0 ), r.NErr( 1 ) )
                              : r.NErr( 0 );
        }

    public:

        static const size_t MAX_N_RES = common::SIntTraits< TLevelCount >::MAX;

        typedef CQueryData TQueryData;

        static bool HaveEqualRanks( const CResult & l, const CResult & r )
        { return (ResultRank( l ) == ResultRank( r )); }

        static bool HaveEqualLevels( const CResult & l, const CResult & r )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return (ResultLevel( l ) == ResultLevel( r ));
        }

        static bool CompareByRank( const CResult & l, const CResult & r )
        { return (ResultRank( l ) > ResultRank( r )); }

        static bool CompareByLevel( const CResult & l, const CResult & r )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return (ResultLevel( l ) < ResultLevel( r ));
        }

        static bool CompareByLevelAtCol( 
                const CResult & l, const CResult & r, size_t col_idx )
        {
            SRPRISM_ASSERT( HaveEqualRanks( l, r ) );
            return (l.NErr( col_idx ) < r.NErr( col_idx ));
        }

        template< bool paired > static int MaxErr( TRank rank )
        { return MAX_N_ERR; }

        static bool HasPairedResult( TRank r )
        { return (r == RANK_PAIRED); }
};

//------------------------------------------------------------------------------
template<> inline bool
CWeakScoringSystem::CQueryData::BestLevelFull< false >( size_t lim ) const
{ return (counts_[0] == lim); }

template<> inline bool
CWeakScoringSystem::CQueryData::BestLevelFull< true >( size_t lim ) const
{
    return (rank_ == RANK_PAIRED && counts_[0] == lim);
}

END_NS( srprism )
END_STD_SCOPES

#endif

