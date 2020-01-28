/*  $Id: out_sam.cpp 536631 2017-05-22 12:51:58Z morgulis $
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
 * File Description: output in SAM format
 *
 */

#include <ncbi_pch.hpp>

#include <vector>
#include <algorithm>

#include "../seq/seqdef.hpp"
#include "out_sam.hpp"

START_STD_SCOPES
START_NS( srprism )

USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
const char * COutSAM::SSAMRecord::STag::VALID_TYPES = "AifzH";

//------------------------------------------------------------------------------
const std::string COutSAM::SEP = "\t";

//------------------------------------------------------------------------------
std::string COutSAM::ComputeMDTag( const CResult & result, int idx )
{
    std::ostringstream os;
    int n_err( result.NErr( idx ) );
    TSeqSize alen( result.GetAlignLen( idx ) ),
             left_offset( result.GetLeftOffset( idx ) );

    if( n_err > 0 ) {
        CResult::CErrorIterator err_i( result.ErrorIterator( idx ) );
        TSeqSize pos_prev( 0 );
        TErrType type_prev( SErrType::E );
        size_t ce_count( 0 ), adj( 0 );

        do {
            TSeqSize pos( 1 + err_i.ErrPos() - left_offset );

            if( err_i.ErrType() != SErrType::I ) {
                if( pos - pos_prev > 1 ) {
                    os << pos - pos_prev - 1 - adj;
                    ce_count = 0;
                    adj = 0;
                }
                else if( err_i.ErrType() != type_prev ) ce_count = 0;

                if( ce_count == 0 && err_i.ErrType() == SErrType::D ) {
                    os << '^';
                }

                TLetter sletter = seq_store_->GetLetter( 
                        result.SNum(), result.SOff( idx ) + err_i.SOff() - left_offset );
                os << sletter;
                ++ce_count;
                pos_prev = pos;
                type_prev = err_i.ErrType();
            }
            else ++adj;

            err_i.Next();
        } while( !err_i.End() );

        if( alen - pos_prev > 0 ) os << alen - pos_prev - adj;
    }
    else os << alen;

    return os.str();
}

//------------------------------------------------------------------------------
namespace {
    struct S_CIGAR_Elmt
    {
        S_CIGAR_Elmt( TSeqSize c = 0, TErrType e = SErrType::E )
            : count( c ), err_type( e )
        {}

        TSeqSize count;
        TErrType err_type;
    };
}

std::string COutSAM::ComputeCIGAR( const CResult & result, int idx )
{
    std::ostringstream os;
    int n_err( result.NErr( idx ) );
    char type_letter[] = { '=', 'X', 'I', 'D' };

    typedef std::vector< S_CIGAR_Elmt > TErrVec;
    TErrVec err_vec;
    S_CIGAR_Elmt curr_elmt;
    Uint2 lo( result.GetLeftOffset( idx ) ),
          ro( result.GetRightOffset( idx ) );
    if( lo > 0 ) os << lo << 'S';

    if( n_err > 0 ) {
        CResult::CErrorIterator err_i( result.ErrorIterator( idx ) );
        TSeqSize pos_prev( 0 );

        do {
            TSeqSize pos( 1 + err_i.ErrPos() );

            if( pos - pos_prev > 1 ) {
                if( pos_prev > 0 ) err_vec.push_back( curr_elmt );

                if( pos_prev > 0 || pos > lo ) {
                    err_vec.push_back( S_CIGAR_Elmt( pos - pos_prev - 1 ) );
                }

                curr_elmt = S_CIGAR_Elmt( 1, err_i.ErrType() );
                pos_prev = pos;
            }
            else {
                if( curr_elmt.err_type == err_i.ErrType() ) {
                    ++curr_elmt.count;
                    ++pos_prev;
                }
                else {
                    if( pos_prev > 0 ) {
                        err_vec.push_back( S_CIGAR_Elmt( curr_elmt ) );
                    }

                    curr_elmt = S_CIGAR_Elmt( 1, err_i.ErrType() );
                    pos_prev = pos;
                }
            }

            err_i.Next();
        } while( !err_i.End() );

        err_vec.push_back( curr_elmt );
        TSeqSize pos( result.GetFullAlignLen( idx ) );

        if( pos - pos_prev > 0 ) {
            err_vec.push_back( S_CIGAR_Elmt( pos - pos_prev ) );
        }

        if( err_vec[0].err_type == SErrType::E ) {
            err_vec[0].count -= lo;
        }

        if( err_vec[err_vec.size() - 1].err_type == SErrType::E ) {
            err_vec[err_vec.size() - 1].count -= ro;
        }

        for( TErrVec::const_iterator i( err_vec.begin() ); 
                i != err_vec.end(); ++i ) {
            if( i->count != 0 ) os << i->count << type_letter[i->err_type];
        }
    }
    else if( result.GetAlignLen( idx ) > 0 ) {
        os << result.GetAlignLen( idx ) << '=';
    }

    if( ro > 0 ) os << ro << 'S';
    return os.str();
}

//------------------------------------------------------------------------------
std::string COutSAM::ComputeSeq( int idx, bool reverse )
{
    std::string result( in_p_->QData( idx ) );

    if( reverse ) {
        typedef SCodingTraits< CODING_IUPACNA > TTraits;
        std::string result_rc;

        for( size_t i( result.size() ); i > 0; ) {
            result_rc.push_back( TTraits::RC[(unsigned char)result[--i]] );
        }

        return result_rc;
    }
    else return result;
}

//------------------------------------------------------------------------------
void COutSAM::SetUpQData( int idx, SSAMRecord & r )
{
    if( input_fmt_ == "cfasta" || input_fmt_ == "cfastq" )
    {
        r.AddTag( "CS", 'Z', IUPACNA2Color( r.seq, 'A' ) );
        r.AddTag( "CQ", 'Z', in_p_->Qual( idx ) );
        r.qstr = "*";
    }
    else r.qstr = in_p_->Qual( idx );
}

//------------------------------------------------------------------------------
void COutSAM::EmptyOut( 
        int idx, Uint4 mpos, const std::string & sname, TStrand mstrand, 
        bool primary )
{
    SSAMRecord r;
    r.extra_tags = extra_tags_;

    if( no_qids_ ) {
        std::ostringstream os;
        os << in_p_->QId() - 1;
        r.qname = os.str();
    }
    else r.qname = in_p_->Id();

    r.flags |= SSAMRecord::SEQ_UNMAPPED_FLAG;
    if( !primary ) r.flags |= SSAMRecord::NOT_PRIMARY_FLAG;

    if( paired_ ) {
        r.flags |= SSAMRecord::PAIRED_QUERY_FLAG;

        if( idx == 0 ) r.flags |= SSAMRecord::SEQ_FIRST_FLAG;
        else           r.flags |= SSAMRecord::SEQ_SECOND_FLAG;
    }

    r.sname = "*";

    if( mpos == 0 ) r.flags |= SSAMRecord::MATE_UNMAPPED_FLAG;
    else {
        if( mstrand == STRAND_RV ) r.flags |= SSAMRecord::MATE_STRAND_FLAG;
        r.msname = sname;
    }

    r.pos = 0;
    r.CIGAR_str = "*";
    r.mpos = mpos;
    r.seq = in_p_->QData( idx );
    r.quality = 0;
    SetUpQData( idx, r );

    if( out_xa_ ) {
        TQNum qn( (TQNum)( (in_p_->QId() - q_adj_ - 1) ) );
        if( paired_ ) qn *= 2;
        qn += idx;
        r.AddITag( "XA", 'i', qs_->HasRepHashes( qn ) ? 0 : 1 );
    }

    (*os_) << r.Format() << std::endl;
}

//------------------------------------------------------------------------------
void COutSAM::ResultOut( 
        const CResult & result, bool mate_unmapped, TQueryOrdId q_adj,
        int const * pg, bool primary )
{
    if( skip_unmapped_ ) in_p_->Skip( 1 + result.QOrdId( q_adj, paired_ ) );
    else {
        TQueryOrdId qid( result.QOrdId( q_adj, paired_ ) + 1 );

        while( in_p_->QId() < qid ) {
            in_p_->Skip( in_p_->QId() + 1 );
            if( in_p_->QId() >= qid ) break;
            EmptyOut( 0 );
            if( paired_ ) EmptyOut( 1 );
        }
    }

    std::string res_qid( in_p_->Id() ), res_sid;

    if( no_qids_ ) {
        std::ostringstream os;
        os << in_p_->QId() - 1; 
        res_qid = os.str();
    }

    if( sid_map_ != 0 ) res_sid = (*sid_map_)[result.SNum()];
    else {
        std::ostringstream os;
        os << result.SNum();
        res_sid = os.str();
    }

    int quality( result.Quality() );

    char mdtag[] = { 'M', 'D' };
    char nmtag[] = { 'N', 'M' };

    // Uint2 lo_1( result.GetLeftOffset( 0 ) );

    if( paired_ ) {
        if( result.Paired() ) {
            // Uint2 lo_2( result.GetLeftOffset( 1 ) );
            SSAMRecord sam_record_1, sam_record_2;
            sam_record_1.extra_tags = sam_record_2.extra_tags = extra_tags_;
            sam_record_1.qname = sam_record_2.qname = res_qid;
            sam_record_1.quality = sam_record_2.quality = quality;
            SSAMRecord::TFlags flags_1( 0 ), flags_2( 0 );
            flags_1 |= SSAMRecord::PAIRED_QUERY_FLAG;
            flags_2 |= SSAMRecord::PAIRED_QUERY_FLAG;
            flags_1 |= SSAMRecord::PAIRED_ALIGN_FLAG;
            flags_2 |= SSAMRecord::PAIRED_ALIGN_FLAG;

            if( !primary ) {
                flags_1 |= SSAMRecord::NOT_PRIMARY_FLAG;
                flags_2 |= SSAMRecord::NOT_PRIMARY_FLAG;
            }

            if( result.Strand( 0 ) == seq::STRAND_RV ) {
                flags_1 |= SSAMRecord::SEQ_STRAND_FLAG;
                flags_2 |= SSAMRecord::MATE_STRAND_FLAG;
            }

            if( result.Strand( 1 ) == seq::STRAND_RV ) {
                flags_2 |= SSAMRecord::SEQ_STRAND_FLAG;
                flags_1 |= SSAMRecord::MATE_STRAND_FLAG;
            }

            Sint8 diff( 0 );

            {
                /*
                Sint8 l1( result.SOff( 0 ) + result.GetLeftOffset( 0 ) ),
                      l2( result.SOff( 1 ) + result.GetLeftOffset( 1 ) );
                l1 -= result.GetLeftOffset( 0 );
                l2 -= result.GetLeftOffset( 1 );
                Sint8 l( std::min( l1, l2 ) ),
                      r( std::max( l1 + result.GetFullSubjLen( 0 ),
                                   l2 + result.GetFullSubjLen( 1 ) ) );
                */
                Sint8 l1( result.SOff( 0 ) ),
                      l2( result.SOff( 1 ) );
                Sint8 l( std::min( l1, l2 ) ),
                      r( std::max( l1 + result.GetSubjLen( 0 ), 
                                   l2 + result.GetSubjLen( 1 ) ) );
                SRPRISM_ASSERT( r >= l );
                diff = r - l;
            }

            flags_1 |= SSAMRecord::SEQ_FIRST_FLAG;
            flags_2 |= SSAMRecord::SEQ_SECOND_FLAG;
            sam_record_1.diff = sam_record_2.diff = diff;
            sam_record_1.flags = flags_1;
            sam_record_2.flags = flags_2;
            sam_record_1.sname = sam_record_2.sname = res_sid;
            sam_record_1.msname = sam_record_2.msname = res_sid;
            sam_record_1.pos = 1 + seq_store_->GetAdjustedPos( 
                    result.SNum(), result.SOff( 0 ) );
                    // result.SNum(), lo_1 + result.SOff( 0 ) );
            sam_record_2.pos = 1 + seq_store_->GetAdjustedPos(
                    result.SNum(), result.SOff( 1 ) );
                    // result.SNum(), lo_2 + result.SOff( 1 ) );
            sam_record_1.CIGAR_str = ComputeCIGAR( result, 0 );
            sam_record_2.CIGAR_str = ComputeCIGAR( result, 1 );
            sam_record_1.mpos = 1 + seq_store_->GetAdjustedPos(
                    result.SNum(), result.SOff( 1 ) );
                    // result.SNum(), lo_2 + result.SOff( 1 ) );
            sam_record_2.mpos = 1 + seq_store_->GetAdjustedPos(
                    result.SNum(), result.SOff( 0 ) );
                    // result.SNum(), lo_1 + result.SOff( 0 ) );
            sam_record_1.seq = ComputeSeq( 
                    0, (result.Strand( 0 ) == seq::STRAND_RV) );
            sam_record_2.seq = ComputeSeq( 
                    1, (result.Strand( 1 ) == seq::STRAND_RV) );
            SetUpQData( 0, sam_record_1 );
            SetUpQData( 1, sam_record_2 );
            sam_record_1.AddTag( mdtag, 'Z', ComputeMDTag( result, 0 ) );
            sam_record_2.AddTag( mdtag, 'Z', ComputeMDTag( result, 1 ) );
            sam_record_1.AddITag( nmtag, 'i', result.NErr( 0 ) );
            sam_record_2.AddITag( nmtag, 'i', result.NErr( 1 ) );

            if( out_xa_ ) {
                sam_record_1.AddITag( "XA", 'i', pg[0] );
                sam_record_2.AddITag( "XA", 'i', pg[1] );
            }

            (*os_) << sam_record_1.Format() << std::endl
                   << sam_record_2.Format() << std::endl;
        }
        else {
            int idx( result.PairPos() );

            if( idx == 1 && mate_unmapped && (primary || !skip_unmapped_) ) {
                EmptyOut( 
                        0, 1 + seq_store_->GetAdjustedPos( 
                            result.SNum(), result.SOff( 0 ) ),
                            // result.SNum(), lo_1 + result.SOff( 0 ) ),
                        res_sid, result.Strand( 0 ), primary );
            }

            SSAMRecord sam_record;
            sam_record.extra_tags = extra_tags_;
            sam_record.qname = res_qid;
            sam_record.quality = quality;
            SSAMRecord::TFlags flags( 0 );
            flags |= SSAMRecord::PAIRED_QUERY_FLAG;
            if( !primary ) flags |= SSAMRecord::NOT_PRIMARY_FLAG;

            if( result.Strand( 0 ) == seq::STRAND_RV ) {
                flags |= SSAMRecord::SEQ_STRAND_FLAG;
            }

            if( mate_unmapped ) flags |= SSAMRecord::MATE_UNMAPPED_FLAG;

            if( idx == 0 ) flags |= SSAMRecord::SEQ_FIRST_FLAG;
            else           flags |= SSAMRecord::SEQ_SECOND_FLAG;

            sam_record.flags = flags;
            sam_record.sname = res_sid;
            sam_record.pos = 1 + seq_store_->GetAdjustedPos( 
                    result.SNum(), result.SOff( 0 ) );
                    // result.SNum(), lo_1 + result.SOff( 0 ) );
            sam_record.CIGAR_str = ComputeCIGAR( result, 0 );
            sam_record.seq = ComputeSeq(
                    idx, (result.Strand( 0 ) == seq::STRAND_RV) );
            SetUpQData( idx, sam_record );
            sam_record.AddTag( mdtag, 'Z', ComputeMDTag( result, 0 ) );
            sam_record.AddITag( nmtag, 'i', result.NErr( 0 ) );
            
            if( out_xa_ ) sam_record.AddITag( "XA", 'i', pg[idx] );

            (*os_) << sam_record.Format() << std::endl;

            if( idx == 0 && mate_unmapped && (primary || !skip_unmapped_) ) {
                EmptyOut( 
                        1, 1 + seq_store_->GetAdjustedPos( 
                            result.SNum(), result.SOff( 0 ) ),
                            // result.SNum(), lo_1 + result.SOff( 0 ) ),
                        res_sid, result.Strand( 0 ), primary );
            }
        }
    }
    else {
        SSAMRecord sam_record;
        sam_record.extra_tags = extra_tags_;
        sam_record.quality = quality;
        sam_record.flags = result.Strand( 0 ) == seq::STRAND_RV 
                         ? SSAMRecord::SEQ_STRAND_FLAG : 0;
        if( !primary ) sam_record.flags |= SSAMRecord::NOT_PRIMARY_FLAG;
        sam_record.qname     = res_qid;
        sam_record.sname     = res_sid;
        sam_record.pos       = 1 + seq_store_->GetAdjustedPos( 
                result.SNum(), result.SOff( 0 ) );
                // result.SNum(), lo_1 + result.SOff( 0 ) );
        sam_record.CIGAR_str = ComputeCIGAR( result, 0 );
        sam_record.seq = 
            ComputeSeq( 0, (result.Strand( 0 ) == seq::STRAND_RV) );
        SetUpQData( 0, sam_record );
        sam_record.AddTag( mdtag, 'Z', ComputeMDTag( result, 0 ) );
        sam_record.AddITag( nmtag, 'i', result.NErr( 0 ) );

        if( out_xa_ ) sam_record.AddITag( "XA", 'i', pg[0] );

        (*os_) << sam_record.Format() << std::endl;
    }
}

//------------------------------------------------------------------------------
void COutSAM::ResultOut( 
        const CResult & result_1, const CResult & result_2, TQueryOrdId q_adj,
        int const * pg )
{
    SRPRISM_ASSERT( paired_ );

    if( skip_unmapped_ ) in_p_->Skip( 1 + result_1.QOrdId( q_adj, paired_ ) );
    else {
        TQueryOrdId qid( result_1.QOrdId( q_adj, paired_ ) + 1 );

        while( in_p_->QId() < qid ) {
            in_p_->Skip( in_p_->QId() + 1 );
            if( in_p_->QId() >= qid ) break;
            EmptyOut( 0 );
            if( paired_ ) EmptyOut( 1 );
        }
    }

    std::string res_qid( in_p_->Id() ), res_sid_1, res_sid_2;

    if( no_qids_ ) {
        std::ostringstream os;
        os << in_p_->QId() - 1;
        res_qid = os.str();
    }

    if( sid_map_ != 0 ) {
        res_sid_1 = (*sid_map_)[result_1.SNum()];
        res_sid_2 = (*sid_map_)[result_2.SNum()];
    }
    else {
        {
            std::ostringstream os;
            os << result_1.SNum();
            res_sid_1 = os.str();
        }
        {
            std::ostringstream os;
            os << result_2.SNum();
            res_sid_2 = os.str();
        }
    }

    int qual_1( result_1.Quality() );
    int qual_2( result_2.Quality() );

    char mdtag[] = { 'M', 'D' };
    char nmtag[] = { 'N', 'M' };

    /*
    Uint2 lo_1( result_1.GetLeftOffset( 0 ) );
    Uint2 lo_2( result_2.GetLeftOffset( 0 ) );
    */
    SSAMRecord sam_record_1, sam_record_2;
    sam_record_1.extra_tags = sam_record_2.extra_tags = extra_tags_;
    sam_record_1.qname = sam_record_2.qname = res_qid;
    sam_record_1.quality = qual_1;
    sam_record_2.quality = qual_2;
    SSAMRecord::TFlags flags_1( 0 ), flags_2( 0 );
    flags_1 |= SSAMRecord::PAIRED_QUERY_FLAG;
    flags_2 |= SSAMRecord::PAIRED_QUERY_FLAG;

    if( result_1.Strand( 0 ) == seq::STRAND_RV ) {
        flags_1 |= SSAMRecord::SEQ_STRAND_FLAG;
        flags_2 |= SSAMRecord::MATE_STRAND_FLAG;
    }

    if( result_2.Strand( 0 ) == seq::STRAND_RV ) {
        flags_2 |= SSAMRecord::SEQ_STRAND_FLAG;
        flags_1 |= SSAMRecord::MATE_STRAND_FLAG;
    }

    flags_1 |= SSAMRecord::SEQ_FIRST_FLAG;
    flags_2 |= SSAMRecord::SEQ_SECOND_FLAG;
    sam_record_1.diff = sam_record_2.diff = 0;
    sam_record_1.flags = flags_1;
    sam_record_2.flags = flags_2;
    sam_record_1.sname = res_sid_1;
    sam_record_2.sname = res_sid_2;
    sam_record_1.msname = res_sid_2;
    sam_record_2.msname = res_sid_1;
    sam_record_1.pos = 1 + seq_store_->GetAdjustedPos( 
            result_1.SNum(), result_1.SOff( 0 ) );
            // result_1.SNum(), lo_1 + result_1.SOff( 0 ) );
    sam_record_2.pos = 1 + seq_store_->GetAdjustedPos(
            result_2.SNum(), result_2.SOff( 0 ) );
            // result_2.SNum(), lo_2 + result_2.SOff( 0 ) );
    sam_record_1.CIGAR_str = ComputeCIGAR( result_1, 0 );
    sam_record_2.CIGAR_str = ComputeCIGAR( result_2, 0 );
    sam_record_1.mpos = sam_record_2.pos;
    sam_record_2.mpos = sam_record_1.pos;
    sam_record_1.seq = ComputeSeq( 
            0, (result_1.Strand( 0 ) == seq::STRAND_RV) );
    sam_record_2.seq = ComputeSeq( 
            1, (result_2.Strand( 0 ) == seq::STRAND_RV) );
    SetUpQData( 0, sam_record_1 );
    SetUpQData( 1, sam_record_2 );
    sam_record_1.AddTag( mdtag, 'Z', ComputeMDTag( result_1, 0 ) );
    sam_record_2.AddTag( mdtag, 'Z', ComputeMDTag( result_2, 0 ) );
    sam_record_1.AddITag( nmtag, 'i', result_1.NErr( 0 ) );
    sam_record_2.AddITag( nmtag, 'i', result_2.NErr( 0 ) );

    if( out_xa_ ) {
        sam_record_1.AddITag( "XA", 'i', pg[0] );
        sam_record_2.AddITag( "XA", 'i', pg[1] );
    }

    (*os_) << sam_record_1.Format() << std::endl
           << sam_record_2.Format() << std::endl;
}

//------------------------------------------------------------------------------
void COutSAM::FinalizeBatch() {
    if( !skip_unmapped_ ) {
        TQueryOrdId last_id( qs_->size() );
        if( paired_ ) last_id /= 2;
        // last_id += q_adj_ + 1;
        last_id += q_adj_;

        while( !in_p_->Done() && in_p_->QId() < last_id ) {
            if( !in_p_->Skip( in_p_->QId() + 1 ) ) break;
            // if( in_p_->QId() >= last_id ) break;
            EmptyOut( 0 );
            if( paired_ ) EmptyOut( 1 );
            // if( in_p_->QId() >= last_id ) break;
        }
    }
}

//------------------------------------------------------------------------------
COutSAM::~COutSAM()
{
    /*
    if( !skip_unmapped_ ) {
        while( !in_p_->Done() ) {
            if( !in_p_->Skip( in_p_->QId() + 1 ) ) break;
            EmptyOut( 0 );
            if( paired_ ) EmptyOut( 1 );
        }
    }
    */
}

END_NS( srprism )
END_STD_SCOPES

