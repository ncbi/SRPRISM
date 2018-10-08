/*  $Id: result.hpp 536631 2017-05-22 12:51:58Z morgulis $
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
 * File Description: representation of a search result
 *
 */

#ifndef __SRPRISM_RESULT_HPP__
#define __SRPRISM_RESULT_HPP__

#include "../common/def.h"

#include <sstream>

#ifndef NCBI_CPP_TK

#include <srprism/srprismdef.hpp>
#include <srprism/seqstore.hpp>
#include <srprism/align.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqstore.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/align.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CResult
{
    private:

        typedef TDBOrdId TSNum;

        struct SResData
        {
            TQNum q_n;
            TSNum s_n;
            common::Sint2 quality;
            common::Uint1 n_col;
        };

        static const size_t RES_DATA_LEN = sizeof( SResData );

        struct SColData
        {
            TSeqSize s_pos;

            struct
            {
                unsigned int left_offset : 14;
                unsigned int right_offset : 14;
                unsigned int strand : 1;
            } f;

            common::Uint2 align_len;
            common::Uint1 n_err;

            common::Uint2 FullAlignLen( void ) const
            { return align_len + f.left_offset + f.right_offset; }
        };

        static const size_t COL_DATA_LEN = sizeof( SColData );

        struct SErrDescr
        {
            unsigned int err_pos : 14;
            unsigned int err_type : 2;
        };

        static const size_t ERR_DESCR_LEN = sizeof( SErrDescr );

        SResData * GetResData( void ) const { return (SResData *)data_; }

        SColData * GetColData( void ) const
        { return (SColData *)(data_ + RES_DATA_LEN); }

        SErrDescr * GetErrDescr( size_t col_idx ) const
        {
            SColData * curr_col( GetColData() );
            SErrDescr * res( 
                    (SErrDescr *)(GetColData() + GetResData()->n_col) );
            while( col_idx-- != 0 ) res += (curr_col++)->n_err;
            return res;
        }

    public:

        class CErrorIterator
        {
            public:

                CErrorIterator( const CResult & r, size_t col_idx )
                    : align_len_( r.GetFullAlignLen( col_idx ) ),
                      s_( r.Strand( col_idx ) ), 
                      adj_( 0 ), iadj_( 0 ), n_del_( r.GetNDel( col_idx ) ),
                      curr_err_( r.GetErrDescr( col_idx ) ),
                      end_err_( curr_err_ + r.GetColData()[col_idx].n_err )
                {
                    if( curr_err_ != end_err_ ) { 
                        if( curr_err_->err_type == SErrType::D ) ++adj_;
                        else if( curr_err_->err_type == SErrType::I ) ++iadj_;
                    }
                }

                bool End( void ) const { return (curr_err_ == end_err_); }

                void Next( void )
                { 
                    ++curr_err_;

                    if( !End() ) {
                        if( curr_err_->err_type == SErrType::D ) ++adj_;
                        else if( curr_err_->err_type == SErrType::I ) ++iadj_;
                    }
                }

                TSeqSize QOff( void ) const { return curr_err_->err_pos; }

                TSeqSize ErrPos( void ) const 
                {
                    TSeqSize qp;

                    if( s_ == seq::STRAND_FW ) qp = curr_err_->err_pos;
                    else { 
                        qp = align_len_ - n_del_ - curr_err_->err_pos - 1;
                        if( ErrType() == SErrType::D ) --qp;
                    }

                    return qp + adj_;
                }

                TSeqSize SOff( void ) const { return ErrPos() - iadj_; }

                TErrType ErrType( void ) const { return curr_err_->err_type; }

            private:

                TSeqSize align_len_;
                TStrand s_;
                common::Uint1 adj_, iadj_, n_del_;
                const SErrDescr * curr_err_, * end_err_;
        };

        static const size_t EstimateLen( 
                size_t n_col, size_t n_err_1, size_t n_err_2 = 0 )
        {
            size_t res( RES_DATA_LEN + n_col*COL_DATA_LEN );
            SRPRISM_ASSERT( res%4 == 0 );
            size_t add( (n_err_1 + n_err_2)*ERR_DESCR_LEN );
            add = (add == 0) ? 0 : ((((add - 1)>>2) + 1)<<2);
            return res + add;
        }

        CErrorIterator ErrorIterator( size_t col_idx ) const
        { return CErrorIterator( *this, col_idx ); }

    private:

        void InitResData( TQNum q_n, TSNum s_n, common::Uint1 n_col )
        {
            GetResData()->q_n = q_n;
            GetResData()->s_n = s_n;
            GetResData()->n_col = n_col;
        }

        void InitColData( 
                size_t col_idx, TSeqSize s_pos, common::Uint2 align_len,
                common::Uint2 left_offset, common::Uint2 right_offset,
                TStrand strand, common::Uint1 n_err )
        {
            GetColData()[col_idx].s_pos = s_pos;
            GetColData()[col_idx].align_len = align_len;
            GetColData()[col_idx].f.left_offset = left_offset;
            GetColData()[col_idx].f.right_offset = right_offset;
            GetColData()[col_idx].f.strand = strand;
            GetColData()[col_idx].n_err = n_err;
        }

    public:

        CResult( char * data ) : data_( data ) {}

        void Init( 
                TQNum q_n, TSNum s_n, 
                TSeqSize s_pos, common::Uint2 align_len, 
                common::Uint2 left_offset, common::Uint2 right_offset,
                TStrand strand, common::Uint1 n_err )
        {
            InitResData( q_n, s_n, 1 );
            InitColData( 
                    0, s_pos, align_len, 
                    left_offset, right_offset, strand, n_err );
        }

        void Init(
                TQNum q_n, TSNum s_n,
                TSeqSize s_pos_1, TSeqSize s_pos_2,
                common::Uint2 align_len_1, common::Uint2 align_len_2,
                common::Uint2 left_offset_1, common::Uint2 left_offset_2,
                common::Uint2 right_offset_1, common::Uint2 right_offset_2,
                TStrand strand_1, TStrand strand_2, 
                common::Uint1 n_err_1, common::Uint1 n_err_2 )
        {
            InitResData( q_n, s_n, 2 );
            InitColData( 
                    0, s_pos_1, align_len_1, 
                    left_offset_1, right_offset_1, strand_1, n_err_1 );
            InitColData( 
                    1, s_pos_2, align_len_2, 
                    left_offset_2, right_offset_2, strand_2, n_err_2 );
        }

        void ReverseStrands( void )
        {
            size_t n_col( GetResData()->n_col );

            for( size_t i( 0 ); i < n_col; ++i ) {
                SColData & sd( GetColData()[i] );
                common::Uint2 fal( sd.FullAlignLen() ),
                              // sal( fal - GetNIns( i ) ),
                              qal( fal - GetNDel( i ) ),
                              al( sd.align_len ),
                              ssal( al - GetNIns( i ) );
                int n_err( sd.n_err );
                SErrDescr * ed( GetErrDescr( i ) );

                // adjust error positions
                for( int i( 0 ); i < n_err; ++i ) {
                    ed[i].err_pos = qal - ed[i].err_pos - 1;
                    if( ed[i].err_type == SErrType::D ) --ed[i].err_pos;
                }

                // flip left and right soft masked regions
                unsigned int t( sd.f.left_offset );
                sd.f.left_offset = sd.f.right_offset;
                sd.f.right_offset = t;

                // reverse the strand
                sd.f.strand = seq::ReverseStrand( sd.f.strand );

                // adjust subject position
                sd.s_pos = sd.f.strand == seq::STRAND_FW ?
                    sd.s_pos - ssal : sd.s_pos + ssal;
                    // sd.s_pos - sal : sd.s_pos + sal;
            }
        }

        void SetErrorInfo( 
                size_t col_idx, size_t err_idx, 
                TSeqSize err_pos, TErrType err_type )
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            SRPRISM_ASSERT( err_idx < GetColData()[col_idx].n_err );
            SErrDescr & err( GetErrDescr( col_idx )[err_idx] );
            err.err_pos = err_pos;
            err.err_type = err_type;
        }

        size_t GetRawLen( void ) const
        {
            size_t n_col( GetResData()->n_col );
            size_t n_err( 0 );

            for( size_t i( 0 ); i < n_col; ++i ) { 
                n_err += GetColData()[i].n_err;
            }

            return EstimateLen( n_col, n_err );
        }

        bool CheckLen( size_t len ) const
        {
            size_t clen( RES_DATA_LEN );
            if( len < clen ) return false;
            size_t n_col( GetResData()->n_col );
            clen += n_col*COL_DATA_LEN;
            if( len < clen ) return false;
            size_t n_err( 0 );

            for( size_t i( 0 ); i < n_col; ++i ) {
                n_err += GetColData()[i].n_err;
            }

            clen = EstimateLen( n_col, n_err );
            return (len >= clen);
        }

        bool Empty( void ) const { return (data_ == 0); }

        void Copy( char * new_loc )
        {
            SRPRISM_ASSERT( new_loc != 0 );
            memcpy( (void *)new_loc, (const void *)data_, GetRawLen() );
            data_ = new_loc;
        }

        void Clone( const CResult & r )
        { memcpy( (void *)data_, (const void *)r.data_, r.GetRawLen() ); }

        TQNum QNum( void ) const { return GetResData()->q_n; }
        void SetQNum( TQNum qn ) { GetResData()->q_n = qn; }

        TDBOrdId SNum( void ) const { return GetResData()->s_n; }

        int GetNId( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            int n_err( GetColData()[col_idx].n_err );
            const SErrDescr * err_descr( GetErrDescr( col_idx ) );
            int res( 0 );

            for( int i( 0 ); i < n_err; ++i, ++err_descr ) {
                if( err_descr->err_type == SErrType::I || 
                        err_descr->err_type == SErrType::D ) {
                    ++res;
                }
            }

            return res;
        }

        int GetNGOpen( size_t col_idx ) const { return 0; }
        
        int GetNIns( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            int n_err( GetColData()[col_idx].n_err );
            const SErrDescr * err_descr( GetErrDescr( col_idx ) );
            int res( 0 );

            for( int i( 0 ); i < n_err; ++i, ++err_descr ) {
                if( err_descr->err_type == SErrType::I ) ++res;
            }

            return res;
        }

        int GetNDel( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            int n_err( GetColData()[col_idx].n_err );
            const SErrDescr * err_descr( GetErrDescr( col_idx ) );
            int res( 0 );

            for( int i( 0 ); i < n_err; ++i, ++err_descr ) {
                if( err_descr->err_type == SErrType::D ) ++res;
            }

            return res;
        }

        int NErr( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].n_err;
        }

        TSeqSize RawSOff( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].s_pos;
        }

        TSeqSize SOff( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            SColData & cd( GetColData()[col_idx] );
            return (cd.f.strand == seq::STRAND_FW) ? 
                cd.s_pos : cd.s_pos + GetNIns( col_idx ) - cd.align_len;
                // cd.s_pos : cd.s_pos + GetNIns( col_idx ) - cd.FullAlignLen();
        }

        bool Paired( void ) const { return (GetResData()->n_col > 1); }

        TQueryOrdId QOrdId( TQueryOrdId adj, bool paired_search ) const
        { return paired_search ? (adj + QNum())/2 : adj + QNum(); }

        void SetQuality( common::Sint2 q ) { GetResData()->quality = q; }
        common::Sint2 Quality( void ) const { return GetResData()->quality; }

        int PairPos(void) const { return (QNum()%2); }

        TStrand Strand( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].f.strand;
        }

        common::Uint2 GetAlignLen( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].align_len;
        }

        common::Uint2 GetSubjLen( size_t col_idx ) const
        {
            return GetAlignLen( col_idx ) - GetNIns( col_idx );
        }

        common::Uint2 GetFullAlignLen( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].FullAlignLen();
        }
        
        common::Uint2 GetFullSubjLen( size_t col_idx ) const
        { return GetFullAlignLen( col_idx ) - GetNIns( col_idx ); }

        common::Uint2 GetRawLeftOffset( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].f.left_offset;
        }

        common::Uint2 GetRawRightOffset( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            return GetColData()[col_idx].f.right_offset;
        }

        common::Uint2 GetLeftOffset( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            SColData & cd( GetColData()[col_idx] );
            return (cd.f.strand == seq::STRAND_FW) ? 
                cd.f.left_offset : cd.f.right_offset;
        }

        common::Uint2 GetRightOffset( size_t col_idx ) const
        {
            SRPRISM_ASSERT( col_idx < GetResData()->n_col );
            SColData & cd( GetColData()[col_idx] );
            return (cd.f.strand == seq::STRAND_FW) ? 
                cd.f.right_offset : cd.f.left_offset;
        }

        struct SHLCompare {
        private:

            bool IsAL( CResult const & r ) const { 
                TDBOrdId sid( r.SNum() );
                return (sid != seqstore->GetRefOId( sid ));
            }

        public:
            SHLCompare( CSeqStore * s ) : seqstore( s ) {}

            bool operator()( CResult const & l, CResult const & r ) const {
                if( l.QNum() == r.QNum() ) {
                    if( IsAL( l ) == IsAL( r ) ) return l.SNum() < r.SNum();
                    else return IsAL( l );
                }
                else return l.QNum() < r.QNum();
            }

        private:

            CSeqStore * seqstore;
        };

        struct SLLCompare {
        private:

            typedef std::pair< TStrand, TStrand > TStrandConf;

            static TStrandConf GetStrandConf( CResult const & r ) {
                if( r.Paired() ) {
                    return std::make_pair( r.Strand( 0 ), r.Strand( 1 ) );
                }
                else return std::make_pair( r.Strand( 0 ), seq::STRAND_FW );
            }

        public:

            static bool IsEquiv( CResult const & l, CResult const & r ) {
                SRPRISM_ASSERT( l.QNum() == r.QNum() );
                SRPRISM_ASSERT( l.Paired() == r.Paired() );
                return GetStrandConf( l ) == GetStrandConf( r ) &&
                    l.SNum() == r.SNum();
            }

            bool operator()( CResult const & l, CResult const & r ) const {
                SRPRISM_ASSERT( l.QNum() == r.QNum() );
                SRPRISM_ASSERT( l.SNum() == r.SNum() );
                SRPRISM_ASSERT( l.Paired() == r.Paired() );

                TStrandConf scl( GetStrandConf( l ) ),
                            scr( GetStrandConf( r ) );

                if( scl == scr ) {
                    if( l.Paired() && l.SOff( 0 ) == r.SOff( 0 ) ) {
                        return l.SOff( 1 ) < r.SOff( 1 );
                    }
                    else return l.SOff( 0 ) < r.SOff( 0 );
                }
                else return scl < scr;
            }
        };

        template< typename t_scoring_sys >
        struct CompareLevels
        {
            typedef t_scoring_sys TScore;

            bool operator()( const CResult & l, const CResult & r ) const
            {
                SRPRISM_ASSERT( l.GetResData()->q_n == r.GetResData()->q_n );
                SRPRISM_ASSERT( TScore::HaveEqualRanks( l, r ) );

                if( TScore::HaveEqualLevels( l, r ) ) {
                    int fine_level_comp( 0 );

                    if( l.Paired() ) {
                        for( size_t i( 0 ); i < 2; ++i ) {
                            if( TScore::CompareByLevelAtCol( l, r, i ) ) {
                                fine_level_comp = -1;
                                break;
                            }

                            if( TScore::CompareByLevelAtCol( r, l, i ) ) {
                                fine_level_comp = 1;
                                break;
                            }
                        }
                    }

                    return (fine_level_comp < 0);
                }
                else return TScore::CompareByLevel( l, r );
            }
        };

    private:

        char * data_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

