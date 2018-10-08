/*  $Id: out_tabular.cpp 431273 2014-04-02 17:10:44Z morgulis $
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
 * File Description: output in tabular format
 *
 */

#include <ncbi_pch.hpp>

#include "../seq/seqdef.hpp"
#include "out_tabular.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( seq )

//------------------------------------------------------------------------------
const std::string COutTabular::SEP = "\t";

//------------------------------------------------------------------------------
void COutTabular::ResultOut( 
        const CResult & result, bool mate_unmapped, TQueryOrdId q_adj, 
        int const * pg, bool primary )
{
    static const TCoding DST_CODING = OUTPUT_CODING;

    in_p_->Skip( 1 + result.QOrdId( q_adj, paired_ ) );
    std::string res_qid( in_p_->Id() ), res_sid;

    if( !paired_ || result.Paired() ) (*os_) << "0";
    else (*os_) << 1 + result.PairPos();

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

    (*os_) << SEP << res_qid << SEP << res_sid;
    int n_aligns = result.Paired() ? 2 : 1;

    for( int i = 0; i < n_aligns; ++i ) {
        int pair_pos( (!paired_ || result.Paired()) ? i : result.PairPos() );
        (*os_) << SEP << 1 + seq_store_->GetAdjustedPos( 
                result.SNum(), result.GetLeftOffset( i ) + result.SOff( i ) );
        (*os_) << SEP << ((result.Strand( i ) == STRAND_FW) ? "+" : "-");
        (*os_) << SEP << result.GetLeftOffset( i )
               << SEP << result.GetRightOffset( i );
        (*os_) << SEP << result.NErr( i );

        if( result.NErr( i ) > 0 ) {
            CResult::CErrorIterator err_i( result.ErrorIterator( i ) );
            TLetter qletter, sletter;

            do{
                switch( err_i.ErrType() ) {
                    case SErrType::E: case SErrType::M:
                        qletter = in_p_->QData( pair_pos )[err_i.QOff()];
                        sletter = seq_store_->GetLetter( 
                                result.SNum(), 
                                result.SOff( i ) + err_i.SOff() );
                        break;

                    case SErrType::I:
                        qletter = in_p_->QData( pair_pos )[err_i.QOff()];
                        sletter = SCodingTraits< DST_CODING >::GAP_LETTER;
                        break;

                    case SErrType::D:
                        qletter = SCodingTraits< DST_CODING >::GAP_LETTER;
                        sletter = seq_store_->GetLetter( 
                                result.SNum(), 
                                result.SOff( i ) + err_i.SOff() );
                        break;

                    default: qletter = sletter = 0; SRPRISM_ASSERT( false );
                }

                if( result.Strand( i ) == STRAND_RV ) {
                    qletter = SCodingTraits< OUTPUT_CODING >::RC[qletter];
                }

                (*os_) << SEP << 1 + err_i.ErrPos() << SEP << qletter << sletter;
                err_i.Next();
            } while( !err_i.End() );
        }
    }

    (*os_) << std::endl;
}

END_NS( srprism )
END_STD_SCOPES

