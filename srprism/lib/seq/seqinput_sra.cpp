/*  $Id: seqinput_sra.cpp 540690 2017-07-10 15:23:27Z morgulis $
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
 */

#include <sstream>

#include <seq/seqinput_sra.hpp>

#ifdef USE_SRA

#include <ncbi-vdb/NGS.hpp>

START_STD_SCOPES
START_NS( seq )

//------------------------------------------------------------------------------
std::string const CSeqInput_SRA::QUAL = "*";

//------------------------------------------------------------------------------
CSeqInput_SRA::CSeqInput_SRA( std::string const & acc, int n_cols )
    : acc_( acc ), n_cols_( n_cols ), n_frags_( 0 ), c_col_( 0 ),
      d0_( s0_, 0 ), d1_( s1_, 0 ),
      run_( ncbi::NGS::openReadCollection( acc_ ) ),
      it_( run_.getReadRange( 1, run_.getReadCount(), Read::all ) ),
      start_( 1 )
{
    if( run_.getReadCount() == 0 )
    {
        done_ = true;
    }
    else
    {
        done_ = false;
    }
}

//------------------------------------------------------------------------------
size_t CSeqInput_SRA::Skip( size_t n )
{
    if( done_ ) return 0;
    start_ += n;
    auto n_reads( run_.getReadCount() );

    if( start_ > n_reads )
    {
        done_ = true;
        return n_reads - (start_ - n);
    }

    it_ = run_.getReadRange( start_, n_reads, Read::all );
    return n;
}

//------------------------------------------------------------------------------
bool CSeqInput_SRA::Next( void )
{
    if( done_ )
    {
        return false;
    }

    c_col_ += n_cols_;

    if( c_col_ < n_frags_ )
    {
        s0_ = s1_;
        d0_.size = d1_.size;
        s1_.clear();
        d1_.size = 0;

        {
            std::ostringstream os;
            os << acc_ << '.' << start_ - 1 << '.' << c_col_ + 1;
            id_ = os.str();
        }

        return true;
    }

    n_frags_ = c_col_ = 0;

    if( !it_.nextRead() )
    {
        done_ = true;
        return false;
    }

    n_frags_ = it_.getNumFragments();

    {
        std::ostringstream os;
        os << acc_ << '.' << start_++;

        if( n_frags_ > n_cols_ )
        {
            os << '.' << c_col_ + 1;
        }

        id_ = os.str();
    }
    
    if( it_.getNumFragments() > 0 )
    {
        it_.nextFragment();
        StringRef sr( it_.getFragmentBases() );
        d0_.size = sr.size();
        s0_.resize( sr.size() );

        for( size_t i( 0 ); i < sr.size(); ++i )
        {
            s0_[i] = sr.data()[i];
        }
    }
    else
    {
        s0_.clear();
        d0_.size = 0;
    }

    if( it_.getNumFragments() > 1 )
    {
        it_.nextFragment();
        StringRef sr( it_.getFragmentBases() );
        d1_.size = sr.size();
        s1_.resize( sr.size() );

        for( size_t i( 0 ); i < sr.size(); ++i )
        {
            s1_[i] = sr.data()[i];
        }
    }
    else
    {
        s1_.clear();
        d1_.size = 0;
    }

    return true;
}

END_NS( seq )
END_STD_SCOPES

#endif

