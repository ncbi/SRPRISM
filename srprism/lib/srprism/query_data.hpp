/*  $Id: query_data.hpp 563221 2018-05-04 14:58:33Z morgulis $
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
 * File Description: unique'ed query data structure
 *
 */

#ifndef __SRPRISM_QUERY_DATA_HPP__
#define __SRPRISM_QUERY_DATA_HPP__

#include <cassert>

#include "../common/def.h"

#ifndef NCBI_CPP_TK

#include <srprism/srprismdef.hpp>
#include <srprism/rmap.hpp>
#include <srprism/seqiter.hpp>

#else

#include <../src/internal/align_toolbox/srprism/lib/srprism/srprismdef.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/rmap.hpp>
#include <../src/internal/align_toolbox/srprism/lib/srprism/seqiter.hpp>

#endif

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CQueryData
{
    public:

        typedef common::Uint2 TQLen;

    private:

        typedef common::Uint4 TAuxWord;
        typedef common::Uint2 TFlags;

        //----------------------------------------------------------------------
        // bit field offsets for aux_data_
        //

        // strand of the prefix
        //
        static const size_t STRAND_BIT    = 0;

        // length of the extension for short queries
        //
        static const size_t EXT_LEN_START = 1;
        static const size_t EXT_LEN_END   = 6;
        static const size_t EXT_LEN_LEN   =
            EXT_LEN_END - EXT_LEN_START;

        // first error type for a 'blown-up' hash
        //
        static const size_t ERR_TYPE_1_START = 6;
        static const size_t ERR_TYPE_1_END   = 8;
        static const size_t ERR_TYPE_1_LEN = 
            ERR_TYPE_1_END - ERR_TYPE_1_START;

        // second error type for a 'blown-up' hash
        //
        static const size_t ERR_TYPE_2_START = 8;
        static const size_t ERR_TYPE_2_END   = 10;
        static const size_t ERR_TYPE_2_LEN = 
            ERR_TYPE_2_END - ERR_TYPE_2_START;

        // position of the first error in a 'blown-up' hash
        //
        static const size_t POS_1_START = 10;
        static const size_t POS_1_END   = 15;
        static const size_t POS_1_LEN   =
            POS_1_END - POS_1_START;

        // position of the second error in a 'blown-up' hash
        //
        static const size_t POS_2_START = 15;
        static const size_t POS_2_END   = 20;
        static const size_t POS_2_LEN   =
            POS_2_END - POS_2_START;

        // extension direction flag in aux_data_:
        //      0 - extension to the left of the hash word;
        //      1 - extension to the right of the hash word
        //
        static const size_t EXT_DIR_BIT = 20;

        // indicator that length of the query is less than 32bp
        //
        static const size_t SHORT_BIT = 21;

        // current hash index
        //
        static const size_t CURR_HASH_IDX_START_BIT = 22;
        static const size_t CURR_HASH_IDX_END_BIT = 26;
        static const size_t CURR_HASH_IDX_N_BITS = 
            CURR_HASH_IDX_END_BIT - CURR_HASH_IDX_START_BIT;

        //----------------------------------------------------------------------
        // bit field offsets for flags_
        //

        // query should not be searched at all
        //
        static const size_t IGNORE_BIT = 0;

        // current hash is palindromic
        //
        static const size_t PALINDROME_BIT = 1;

        // query contains ambiguities
        //
        static const size_t AMBIG_BIT = 2;

        //----------------------------------------------------------------------
        // bit field offsets for hash configuration word
        //

        // offset of the start of the seeding area (13 bits)
        //
        static const size_t SA_OFF_START_BIT = 0;
        static const size_t SA_OFF_N_BITS = 
            common::SBinLog< MAX_QUERY_LEN >::VALUE;
        static const size_t SA_OFF_END_BIT = SA_OFF_START_BIT + SA_OFF_N_BITS;

        // number of hashes in the seeding area - 1 (4 bits)
        //
        static const size_t SA_N_HASHES_START_BIT = SA_OFF_END_BIT;
        static const size_t SA_N_HASHES_N_BITS =
            common::SBinLog< MAX_N_HASHES >::VALUE;
        static const size_t SA_N_HASHES_END_BIT = 
            SA_N_HASHES_START_BIT + SA_N_HASHES_N_BITS;

        // for short queries: offset of hash 1;
        // for long queries: offset of hash 2 - HASH_LEN (5 bits)
        //
        static const size_t SA_LAST_HASH_OFF_START_BIT = SA_N_HASHES_END_BIT;
        static const size_t SA_LAST_HASH_OFF_N_BITS = 
            common::SBinLog< HASH_LEN >::VALUE + 1;
        static const size_t SA_LAST_HASH_OFF_END_BIT = 
            SA_LAST_HASH_OFF_START_BIT + SA_LAST_HASH_OFF_N_BITS;

        // bit indicating long SA (overlaps with the previous field)
        //
        static const size_t SA_LONG_BIT = SA_LAST_HASH_OFF_END_BIT - 1;

        // index of the blow up hash for medium and short queries (2 bits)
        //
        static const size_t SA_BLOWUP_HASH_IDX_START_BIT = 
            SA_LAST_HASH_OFF_END_BIT;
        static const size_t SA_BLOWUP_HASH_IDX_N_BITS = 2;
        static const size_t SA_BLOWUP_HASH_IDX_END_BIT =
            SA_BLOWUP_HASH_IDX_START_BIT + SA_BLOWUP_HASH_IDX_N_BITS;

        // bit indicating whether letters should be borrowed from 
        // the right when 'blowup' of a hash occurs
        //
        static const size_t SA_BLOWUP_DIR_BIT = 
            SA_BLOWUP_HASH_IDX_END_BIT;

        // 1 iff SA has ambiguities
        //
        static const size_t SA_AMBIG_BIT = SA_BLOWUP_DIR_BIT + 1;

        //----------------------------------------------------------------------
        // bit field offsets for hash masks
        //

        // bad hash extension direction bit mask: for each hash
        // selects the direction of the extension (1 for right)
        // in the case hash is a 'bad' 16-mer
        //
        static const size_t SA_EXT_DIR_MASK_START_BIT = 0;
        static const size_t SA_EXT_DIR_MASK_END_BIT = 
            SA_EXT_DIR_MASK_START_BIT + MAX_N_HASHES;

        // used bad hash mask: set to 1 for each bad hash that has
        // already been searched.
        //
        static const size_t SA_BAD_MASK_START_BIT = SA_EXT_DIR_MASK_END_BIT;
        static const size_t SA_BAD_MASK_END_BIT = 
            SA_BAD_MASK_START_BIT + MAX_N_HASHES;

        //----------------------------------------------------------------------
        class CBlowUpIterator_Base
        {
            protected:

                static const seq::TCoding CODING = SEQDATA_CODING;
                static const seq::TAlphabet & alphabet;

            public:

                CBlowUpIterator_Base(
                        const CQueryData & q,
                        TSeqSize hash_start, TSeqSize start, TSeqSize end,
                        bool right_ext );

                virtual ~CBlowUpIterator_Base() {}

                virtual bool Next( void ) = 0;
                bool End( void ) const { return done_; }

                TSeqSize GetPos_1( void ) const
                {
                    if( two_err_ && !right_ext_ ) return 0;
                    else return curr_pos_ + curr_pos_adj_;
                }

                TSeqSize GetPos_2( void ) const
                {
                    return curr_pos_adj_ + 
                        ((two_err_ && right_ext_) ? HASH_LEN : curr_pos_);
                }

                TErrType GetType_1( void ) const { return type_1_; }
                TErrType GetType_2( void ) const { return type_2_; }

                TWord GetHashValue( void ) const { return hash_; }

            protected:

                const CQueryData * q_;
                TSeqSize hash_start_, start_, end_;
                TWord hash_, orig_hash_;
                TErrType type_1_, type_2_;
                const bool right_ext_;
                bool done_;
                bool two_err_;
                TSeqSize curr_pos_;
                int curr_pos_adj_;
        };

        class CBlowUpIterator_M : public CBlowUpIterator_Base
        {
            public:

                CBlowUpIterator_M(
                        const CQueryData & q,
                        TSeqSize hash_start, TSeqSize start, TSeqSize end,
                        bool right_ext, bool ambig );

                virtual bool Next( void );

            private:

                seq::TAlphabet::const_iterator alpha_;
                TLetter qletter_;
                bool ambig_;
        };

        class CBlowUpIterator_I : public CBlowUpIterator_Base
        {
            public:

                CBlowUpIterator_I(
                        const CQueryData & q,
                        TSeqSize hash_start, TSeqSize start, TSeqSize end,
                        bool right_ext, bool ambig );

                virtual bool Next( void );

            private:

                static const size_t LETTER_BITS = 
                    seq::SCodingTraits< CODING >::LETTER_BITS;

                static const int ONE_ERR = 0;
                static const int TWO_ERR = 1;

                void RemoveAndShiftLeft( void );
                void RemoveAndShiftRight( void );

                void SearchRight( void );
                void SearchLeft( void );

                seq::TAlphabet::const_iterator alpha_;
                TLetter prev_letter_;
                TLetter qletter_1_, qletter_2_;
        };

        class CBlowUpIterator_D : public CBlowUpIterator_Base
        {
            public:

                CBlowUpIterator_D(
                        const CQueryData & q,
                        TSeqSize hash_start, TSeqSize start, TSeqSize end,
                        bool right_ext, bool ambig );

                virtual bool Next( void );

            private:

                static const size_t LETTER_BITS = 
                    seq::SCodingTraits< CODING >::LETTER_BITS;

                void InsertAndShiftRight( void );
                void InsertAndShiftLeft( void );

                seq::TAlphabet::const_iterator alpha_;
                TLetter prev_letter_, curr_letter_;
        };

        class CBlowUpIterator
        {
            public:

                CBlowUpIterator(
                        const CQueryData & q,
                        TSeqSize hash_start, TSeqSize start, TSeqSize end,
                        bool right_ext, bool ambig );

                bool Next( void );

                bool End( void ) const
                { return curr_iter_ == &iter_d_ && curr_iter_->End(); }

                TSeqSize GetPos_1( void ) const 
                { return curr_iter_->GetPos_1(); }

                TSeqSize GetPos_2( void ) const
                { return curr_iter_->GetPos_2(); }

                TErrType GetType_1( void ) const 
                { return curr_iter_->GetType_1(); }

                TErrType GetType_2( void ) const 
                { return curr_iter_->GetType_2(); }

                TWord GetHashValue( void ) const
                { return curr_iter_->GetHashValue(); }

            private:

                CBlowUpIterator_M iter_m_;
                CBlowUpIterator_I iter_i_;
                CBlowUpIterator_D iter_d_;
                CBlowUpIterator_Base * curr_iter_;
        };

    //--------------------------------------------------------------------------

    public:

        void SetExtLen()
        {
            common::SetField_Unsafe< EXT_LEN_START, EXT_LEN_END >(
                    aux_data_, ExtLen() );
        }

        typedef CBlowUpIterator TBlowUpIterator;

        TBlowUpIterator GetBlowUpIterator( 
                TSeqSize hash_start, TSeqSize start, TSeqSize end, 
                bool right_ext, bool ambig ) const
        { 
            return CBlowUpIterator( 
                    *this, hash_start, start, end, right_ext, ambig ); 
        }

        TSeqSize GetOff( void ) const
        {
            int hash_idx( GetCurrHashIdx() );
            bool BU_hash( !IsSALong() && hash_idx == 3 );
            if( BU_hash ) hash_idx = GetBUHashIdx();
            TSeqSize off( GetHashOffset( hash_idx ) );
            if( BU_hash && !RightBUDir() ) off += NDel() - NIns();
            return off;
        }

        TSeqSize GetExtOff( void ) const
        {
            int hash_idx( GetCurrHashIdx() );
            TSeqSize hash_len( HASH_LEN );

            if( !IsSALong() && hash_idx == 3 ) {
                hash_idx = GetBUHashIdx();

                if( IsRightExtDir() == RightBUDir() ) {
                    hash_len = HASH_LEN + NIns() - NDel();
                }
            }

            TSeqSize hash_start( GetHashOffset( hash_idx ) );
            return IsRightExtDir() 
                        ? hash_start + hash_len
                        : hash_start - std::min( hash_start, hash_len );
        }

        TWord GetExtension( void ) const
        {
            TSeqSize off( GetExtOff() ), ext_len( ExtLen() );
            CSeqFwIterator< SEQDATA_CODING, TWord, TWord > seqiter(
                    Data(), off, off + ext_len );
            TWord ext( (*seqiter).first );

            if( !IsRightExtDir() ) {
                seq::ReverseComplement< SEQDATA_CODING >( ext, ext, ext_len );
            }
            
            return ext;
        }

        class CCompare
        {
            public:

                bool operator()( const CQueryData & l, const CQueryData & r )
                {
                    if( l.Ignored() )
                        if( r.Ignored() ) return (l.QNum() < r.QNum());
                        else return false;
                    else if( r.Ignored() ) return true;
                    else if( l.Prefix() == r.Prefix() ) {
                        if( l.aux_data_ == r.aux_data_ ) {
                            return l.GetExtension() < r.GetExtension();
                            /*
                            if( l.Len() == r.Len() ) {
                                if( l.Data() == r.Data() ) {
                                    return (l.GetCurrHashIdx() < 
                                                r.GetCurrHashIdx());
                                }
                                else return (l.Data() < r.Data());
                            }
                            else return (l.Len() < r.Len());
                            */
                        }
                        else return l.aux_data_ < r.aux_data_;
                    }
                    else return (l.Prefix() < r.Prefix());
                }
        };

        class CCompareRaw
        {
            private:

                typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
                static const size_t WLETTERS = 
                    sizeof( TWord )*TTraits::PACK_FACTOR;
                static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;

            public:

                bool operator()( const CQueryData & l, const CQueryData & r )
                {
                    if( l.Ignored() ) {
                        if( r.Ignored() ) return (l.QNum() < r.QNum());
                        else return false;
                    }
                    else if( r.Ignored() ) return true;
                    else if( l.Len() == r.Len() ) {
                        size_t n_words( 3 + ((l.Len() - 1)>>WSHIFT) );
                        n_words <<= 1;
    
                        for( size_t i = 1; i < n_words - 1; ++i ) {
                            if( l.raw_data_[i] < r.raw_data_[i] ) {
                                return true;
                            }

                            if( l.raw_data_[i] > r.raw_data_[i] ) {
                                return false;
                            }
                        }

                        return (l.QNum() < r.QNum());
                    }
                    else return (l.Len() < r.Len());
                }
        };

        class CEqualRaw
        {
            private:

                typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
                static const size_t WLETTERS = 
                    sizeof( TWord )*TTraits::PACK_FACTOR;
                static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;

            public:

                bool operator()( const CQueryData & l, const CQueryData & r )
                {
                    if( l.Len() == r.Len() ) {
                        size_t n_words( 3 + ((l.Len() - 1)>>WSHIFT) );
                        n_words <<= 1;
    
                        for( size_t i = 1; i < n_words - 1; ++i ) {
                            if( l.raw_data_[i] != r.raw_data_[i] ) return false;
                        }

                        return true;
                    }

                    return false;
                }
        };

        static bool Equiv( const CQueryData & l, const CQueryData & r )
        { return l.Prefix() == r.Prefix() && l.aux_data_ == r.aux_data_; }

        CQueryData( TQNum qnum, TSeqSize len, TWord * raw_data )
            : raw_data_( raw_data ), prefix_( 0 ), aux_data_( 0 ), 
              qnum_( qnum ), len_( 0 ), flags_( 0 )
        { SetLen( len ); }

        TSeqSize FirstHashAmbigPos( TSeqSize qoff ) const
        {
            typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
            static const size_t LBITS = TTraits::LETTER_BITS;
            static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

            SRPRISM_ASSERT( qoff + HASH_LEN <= Len() );
            TWord w( seq::GetWord< SEQDATA_CODING >( Mask(), qoff ) );
            return (common::FirstSetBit_Left( w )>>LSHIFT);
        }

        int NHashAmbigs( TSeqSize qoff ) const
        {
            typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
            static const size_t LBITS = TTraits::LETTER_BITS;
            static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

            SRPRISM_ASSERT( qoff < Len() );
            TWord w( seq::GetWord< SEQDATA_CODING >( Mask(), qoff ) );
            w &= (TWord)0x55555555;
            TSeqSize l( Len() - qoff );
            if( l < HASH_LEN ) w >>= ((HASH_LEN - l)<<LSHIFT);
            return common::CountBits( w );
        }

        TWord GetAmbigMask( TSeqSize qoff ) const
        {
            SRPRISM_ASSERT( qoff < Len() );
            TWord w( seq::GetWord< SEQDATA_CODING >( Mask(), qoff ) );
            return (w&(TWord)0xaaaaaaaa);
        }

        int NAmbigs( TSeqSize start, TSeqSize end ) const
        {
            typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
            static const size_t LBITS = TTraits::LETTER_BITS;
            static const size_t LSHIFT = common::SBinLog< LBITS >::VALUE;

            SRPRISM_ASSERT( start < end );
            SRPRISM_ASSERT( end <= Len() );
            const TWord * m( Mask() );
            int res( 0 );

            while( start + HASH_LEN < end ) {
                TWord w( seq::GetWord< SEQDATA_CODING >( m, start ) );
                w &= (TWord)0x55555555;
                res += common::CountBits( w );
                start += HASH_LEN;
            }

            if( start < end ) {
                TWord w( seq::GetWord< SEQDATA_CODING >( m, start ) );
                TWord mask( common::MaskFront< TWord >( 
                            (end - start)<<LSHIFT ) );
                w &= (mask&0x55555555);
                res += common::CountBits( w );
            }

            return res;
        }

        TPrefix Prefix( void ) const { return prefix_; }
        TQNum QNum( void ) const { return qnum_; }
        void SetQNum( TQNum qn ) { qnum_ = qn; }
        TSeqSize Len( void ) const { return (TSeqSize)len_; }

        void SetPrefix( TPrefix prefix ) { prefix_ = prefix; }

        void SetLen( TSeqSize len ) 
        { 
            SRPRISM_ASSERT( len <= (TSeqSize)common::SIntTraits< TQLen >::MAX );
            len_ = len; 
            common::AssignBit< SHORT_BIT >( 
                    aux_data_, (len_ < MIN_MED_QUERY_LEN) );
        }

        TStrand Strand( void ) const
        { return (TStrand)common::GetBit< STRAND_BIT >( aux_data_ ); }

        void SetStrand( TStrand s ) 
		{ 
			common::AssignBit< STRAND_BIT >( 
				aux_data_, (s == seq::STRAND_RV) ); 
		}

        bool IsRightExtDir( void ) const
        { return common::GetBit< EXT_DIR_BIT >( aux_data_ ); }

        void SetRightExtDir( bool right_ext )
        { common::AssignBit< EXT_DIR_BIT >( aux_data_, right_ext ); }

        bool IsShort( void ) const
        { return (Len() < MIN_MED_QUERY_LEN); }

        TSeqSize ExtLen( void ) const
        { 
            if( IsShort() && IsSALong() ) {
                if( IsRightExtDir() ) return Len() - GetHashOffset( 0 ) - HASH_LEN;
                else return GetHashOffset( 0 );
            }
            else return std::min( EXT_LEN, Len() - HASH_LEN + NDel() - NIns() );
        }

        int GetPos1( void ) const
        {
            return (int)common::GetField< 
                POS_1_START, POS_1_END >( aux_data_ ) - 1;
        }

        void SetPos1( TSeqSize pos )
        {
            SRPRISM_ASSERT( (pos <= 1U + common::SBitFieldTraits< 
                        TSeqSize, POS_1_LEN >::MAX) );
            common::SetField_Unsafe< POS_1_START, POS_1_END >( 
                    aux_data_, (TAuxWord)pos );
        }

        int GetPos2( void ) 
        { 
            return (int)common::GetField< 
                POS_2_START, POS_2_END >( aux_data_ ) - 1;
        }

        void SetPos2( TSeqSize pos )
        {
            SRPRISM_ASSERT( (pos <= 1U + common::SBitFieldTraits< 
                        TSeqSize, POS_2_LEN >::MAX) );
            common::SetField_Unsafe< POS_2_START, POS_2_END >( 
                    aux_data_, (TAuxWord)pos );
        }

        bool Ignored( void ) const
        { return common::GetBit< IGNORE_BIT >( flags_ ); }

        void SetIgnored( bool ignore )
        { common::AssignBit< IGNORE_BIT >( flags_, ignore ); }

        bool IsPalindrome( void ) const
        { return common::GetBit< PALINDROME_BIT >( flags_ ); }

        void SetPalindrome( bool palindrome )
        { common::AssignBit< PALINDROME_BIT >( flags_, palindrome ); }

        bool IsAmbig( void ) const
        { return common::GetBit< AMBIG_BIT >( flags_ ); }

        bool RightExtDir( int idx ) const
        { 
            return common::GetBit( 
                    SA_EXT_DIR_MASK_START_BIT + idx, *(raw_data_ + 1) );
        }

        TSeqSize GetAmbigBUPos( void ) const
        { return FirstHashAmbigPos( GetHashOffset( GetBUHashIdx() ) ); }

        int GetNHashes( void ) const
        { 
            return 1 + common::GetField< 
                SA_N_HASHES_START_BIT, SA_N_HASHES_END_BIT >( *raw_data_ );
        }

        bool HashAmbig( int hash_idx ) const
        { 
            TSeqSize qoff( GetHashOffset( hash_idx ) );
            TWord w( seq::GetWord< SEQDATA_CODING >( Mask(), qoff ) );
            return (w != 0);
        }

        int GetBUHashIdx( void ) const
        { 
            return common::GetField< 
                SA_BLOWUP_HASH_IDX_START_BIT, SA_BLOWUP_HASH_IDX_END_BIT >(
                        *raw_data_ );
        }

        bool IsSALong( void ) const
        { return common::GetBit< SA_LONG_BIT >( *raw_data_ ); }

        int GetCurrHashIdx( void ) const
        { 
            return (int)common::GetField< 
                CURR_HASH_IDX_START_BIT, CURR_HASH_IDX_END_BIT >( aux_data_ );
        }

        bool AmbigSA( void ) const 
        { return common::GetBit< SA_AMBIG_BIT >( *raw_data_ ); }

        TSeqSize GetSAOffset( void ) const
        {
            return common::GetField<
                SA_OFF_START_BIT, SA_OFF_END_BIT >( *raw_data_ );
        }

        TSeqSize GetHashOffset( int hash_idx ) const
        {
            int n_hashes( GetNHashes() );
            bool is_long( IsSALong() );

            if( !is_long && n_hashes < 4 && hash_idx == 3 ) {
                hash_idx = GetBUHashIdx();
            }

            SRPRISM_ASSERT( hash_idx < n_hashes );

            if( hash_idx == 0 || is_long ) {
                return GetSAOffset() + hash_idx*HASH_LEN;
            }
            else if( n_hashes < 3 ) {
                return GetSAOffset() + common::GetField< 
                    SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >(
                            *raw_data_ );
            }
            else if( hash_idx == 1 ) {
                if( RightExtDir( 1 ) ) {
                    return GetSAOffset() + common::GetField< 
                        SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >( 
                                *raw_data_ );
                }
                else return GetSAOffset() + HASH_LEN;
            }
            else {
                return GetSAOffset() + common::GetField< 
                    SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >(
                            *raw_data_ ) + HASH_LEN;
            }
        }

        bool RightExtDir( void ) const
        { 
            int hash_idx( GetCurrHashIdx() );

            if( !IsSALong() && GetNHashes() < 4 && hash_idx == 3 ) {
                hash_idx = GetBUHashIdx();
            }

            return RightExtDir( hash_idx );
        }

    public:

        TSeqSize GetSALen( void ) const
        {
            if( IsSALong() ) return GetNHashes()*HASH_LEN;
            else {
                TSeqSize last_hash_end( common::GetField< 
                        SA_LAST_HASH_OFF_START_BIT, SA_LAST_HASH_OFF_END_BIT >( 
                            *raw_data_ ) + HASH_LEN );
                if( GetNHashes() > 2 ) last_hash_end += HASH_LEN;
                return last_hash_end;
            }
        }

        bool RightBUDir( void ) const
        { return common::GetBit< SA_BLOWUP_DIR_BIT >( *raw_data_ ); }

        bool GetBlowUpParams( TSeqSize & start, TSeqSize & end ) const
        {
            if( IsSALong() ) return false;
            if( GetNHashes() < 2 ) { start = 0; end = HASH_LEN; return true; }
            int bu_idx( GetBUHashIdx() );
            TSeqSize h0( GetHashOffset( 0 ) ), h1( GetHashOffset( 1 ) );

            if( GetSALen() < MIN_MED_QUERY_LEN ) {
                if( bu_idx == 0 ) { start = h1 - h0; end = HASH_LEN; }
                else { start = 0; end = h0 + HASH_LEN - h1; }
            }
            else {
                TSeqSize h2( GetHashOffset( 2 ) );

                switch( bu_idx ) {
                    case 0: start = h1 - h0; end = HASH_LEN; break;
                    case 1: 

                        if( h1 - h0 == HASH_LEN ) {
                            start = h2 - h1;
                            end = HASH_LEN;
                        }
                        else {
                            start = 0;
                            end = h0 + HASH_LEN - h1;
                        }

                        break;

                    case 2: start = 0; end = h1 + HASH_LEN - h2; break;
                    default: SRPRISM_ASSERT( false );
                }
            }

            return true;
        }

        void SetCurrHashIdx( int hash_idx )
        {
            common::SetField< 
                CURR_HASH_IDX_START_BIT, CURR_HASH_IDX_END_BIT >( 
                        aux_data_, (TAuxWord)hash_idx );
        }

        void SetAmbig( bool ambig )
        { common::AssignBit< AMBIG_BIT >( flags_, ambig ); }

        bool GetBadHash( int hash_idx ) const
        {
            return common::GetBit( 
                    SA_BAD_MASK_START_BIT + hash_idx, *(raw_data_ + 1) );
        }

        void SetBadHash( int hash_idx, bool bad )
        { 
            common::AssignBit( 
                    SA_BAD_MASK_START_BIT + hash_idx, 
                    *(raw_data_ + 1), bad );
        }

        void SetBadHash( void )
        { SetBadHash( GetCurrHashIdx(), true ); }

    private:

        void SetSALong( bool sa_long )
        { common::AssignBit< SA_LONG_BIT >( *raw_data_, sa_long ); }

    public:

        const TWord * RawData( void ) const { return raw_data_; }
        void SetRawData( TWord * raw_data ) { raw_data_ = raw_data; }
        const TWord * Data( void ) const { return raw_data_ + 2; }

        const TWord * Mask( void ) const
        {
            typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
            static const size_t WLETTERS = sizeof( TWord )*TTraits::PACK_FACTOR;
            static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;

            size_t n_words( 1 + ((Len() - 1)>>WSHIFT) );
            return Data() + n_words + 1;
        }

    private:

        static common::Uint8 GetWordBadness( const CRMap & rmap, TWord w )
        {
            TWord rw( 0 );
            seq::ReverseComplement< SEQDATA_CODING >( rw, w );
            w = std::min( w, rw );
            return (common::Uint8)(1ULL<<rmap.RepeatRank( w ));
        }

        bool ComputeHC_NormalFixedLen( 
                const CRMap & rmap,
                TSeqSize sa_start, TSeqSize sa_end, 
                TSeqSize start, TSeqSize end, 
                common::Uint8 min_badness,
                const common::Uint8 * bvec,
                const common::Uint8 * mvec,
                const common::Uint8 * avec,
                TWord & hc, TWord & mw, common::Uint8 & badness );

        bool ComputeHC_Normal( 
                const CRMap & rmap,
                TSeqSize sa_start, TSeqSize sa_end, int n_err );

        bool ComputeHC_Short( 
                const CRMap & rmap,
                TSeqSize sa_start, TSeqSize sa_end );

        bool ComputeHC_Med( 
                const CRMap & rmap,
                TSeqSize sa_start, TSeqSize sa_end );

    public:

        bool ComputeHC( 
                const CRMap & rmap,
                TSeqSize sa_start, TSeqSize sa_end, int n_err,
                bool use_fixed_hc, TWord fixed_hc )
        {
            if( use_fixed_hc ) *raw_data_ = fixed_hc;
            else {
                sa_end = std::max( sa_end, MIN_MED_QUERY_LEN );
                sa_end = std::min( sa_end, Len() );

                if( sa_end >= MIN_MED_QUERY_LEN )
                {
                    sa_start = std::min( sa_start, sa_end - MIN_MED_QUERY_LEN );
                }
                else
                {
                    sa_start = 0;
                }

                if( sa_end - sa_start >= (n_err + 1)*HASH_LEN ) {
                    return ComputeHC_Normal( rmap, sa_start, sa_end, n_err );
                }
                else if( n_err == 1 ) {
                    return ComputeHC_Short( rmap, sa_start, sa_end );
                }
                else if( n_err == 2 ) {
                    return ComputeHC_Med( rmap, sa_start, sa_end );
                }
                else SRPRISM_ASSERT( false );

                return false;
            }

            return true;
        }

        int GetSeedNErr( void ) const
        { 
            if( IsShort() ) return MAX_SHORT_ERR;
            else if( !IsSALong() ) return MAX_MED_SEED_N_ERR;
            else return GetNHashes() - 1;
        }

        static size_t NWords( TQLen len )
        {
            typedef seq::SCodingTraits< SEQDATA_CODING > TTraits;
            static const size_t WLETTERS = sizeof( TWord )*TTraits::PACK_FACTOR;
            static const size_t WSHIFT = common::SBinLog< WLETTERS >::VALUE;
            return 1 + ((len - 1)>>WSHIFT);
        }

        size_t NTriplets( void ) const;
        size_t NWords( void ) const { return NWords( Len() ); }

        void SetType1( TErrType t )
        { 
            common::SetField< ERR_TYPE_1_START, ERR_TYPE_1_END >( 
                    aux_data_, (TAuxWord)t );
        }

        TErrType Type1( void ) const
        { 
            return common::GetField< ERR_TYPE_1_START, ERR_TYPE_1_END >( 
                    aux_data_ );
        }

        void SetType2( TErrType t )
        { 
            common::SetField< ERR_TYPE_2_START, ERR_TYPE_2_END >( 
                    aux_data_, (TAuxWord)t );
        }

        TErrType Type2( void ) const
        { 
            return common::GetField< ERR_TYPE_2_START, ERR_TYPE_2_END >( 
                    aux_data_ );
        }

        int NErr( void ) const
        {
            int res( 0 );
            if( Type1() != SErrType::E ) ++res;
            if( Type2() != SErrType::E ) ++res;
            return res;
        }

        int NDel( void ) const
        {
            int res( 0 );
            if( Type1() == SErrType::D ) ++res;
            if( Type2() == SErrType::D ) ++res;
            return res;
        }

        int NIns( void ) const
        {
            int res( 0 );
            if( Type1() == SErrType::I ) ++res;
            if( Type2() == SErrType::I ) ++res;
            return res;
        }

    private:

        TWord * raw_data_;
        TPrefix prefix_;
        TAuxWord aux_data_;
        TQNum qnum_;
        TQLen len_;
        TFlags flags_;
};

END_NS( srprism )
END_STD_SCOPES

#endif

