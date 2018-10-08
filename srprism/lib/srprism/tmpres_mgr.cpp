/*  $Id: tmpres_mgr.cpp 358396 2012-04-02 16:27:41Z morgulis $
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
 * File Description: temporary result storage manipulation
 *
 */

#include <ncbi_pch.hpp>

#include "../common/def.h"

#include "tmpres_mgr.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )

//------------------------------------------------------------------------------
CTmpResMgr::CTmpResBuf::CTmpResBuf( void * buf, size_t buf_size )
    : buf_( (char *)buf ), bufsize_( buf_size ), curr_size_( 0 ), read_pos_( 0 )
{
    std::fill( buf_, buf_ + bufsize_, 0 );
}

//------------------------------------------------------------------------------
void CTmpResMgr::CTmpResBuf::Dump( CWriteBinFile & os )
{
    os.Write( (const char *)buf_, curr_size_ );
    curr_size_ = 0;
}

//------------------------------------------------------------------------------
bool CTmpResMgr::CTmpResBuf::Last(void) const
{
    CResult r( (char *)buf_ + read_pos_ );
    return !r.CheckLen( curr_size_ - read_pos_ );
}

//------------------------------------------------------------------------------
bool CTmpResMgr::CTmpResBuf::Load( common::CReadBinFile & is )
{
    memcpy( 
            (void *)buf_, 
            (const void *)(buf_ + read_pos_), 
            curr_size_ - read_pos_ );
    curr_size_ -= read_pos_; read_pos_ = 0;
    CReadBinFile::TSize bytes( 
            is.Read( (char *)(buf_ + curr_size_), bufsize_ - curr_size_ ) );
    curr_size_ += bytes;
    return (curr_size_ != 0);
}

//------------------------------------------------------------------------------
CResult CTmpResMgr::CTmpResBuf::ReadNext( void )
{ 
    CResult res( (char *)buf_ + read_pos_ ); 
    read_pos_ += res.GetRawLen();
    return res;
}

//------------------------------------------------------------------------------
CTmpResMgr::CTmpResMgr( 
        void * mainbuf, size_t mainbuf_size, const std::string & tmp_name, 
        CTmpStore & tmp_store )
    : tmp_store_( tmp_store ), mainbuf_( mainbuf, mainbuf_size ), 
      tmp_name_( tmp_name ), os_( 0 ), is_( 0 )
{
}

//------------------------------------------------------------------------------
void CTmpResMgr::CheckWriteTmpFile(void)
{
    if( os_.get() == 0 ) {
        os_.reset( new CWriteBinFile( tmp_store_.Register( tmp_name_ ) ) );
    }
}

//------------------------------------------------------------------------------
bool CTmpResMgr::CheckReadTmpFile(void)
{
    if( is_.get() == 0 ) {
        if( !tmp_store_.Find( tmp_name_ ) ) return false;
        is_.reset( new CReadBinFile( tmp_store_.Register( tmp_name_ ) ) );
    }
    
    return true;
}

//------------------------------------------------------------------------------
void CTmpResMgr::LoadInit(void)
{
    if( os_.get() != 0 ) mainbuf_.Dump( *os_ );
    os_.reset( 0 );
    is_.reset( 0 );
    mainbuf_.ReadInit();
}

//------------------------------------------------------------------------------
void CTmpResMgr::LoadFinal(void)
{
    is_.reset( 0 );
    mainbuf_.WriteInit();
}

//------------------------------------------------------------------------------
CResult CTmpResMgr::Load( void )
{
    if( mainbuf_.Last() ) {
        if( !CheckReadTmpFile() ) return CResult( 0 );
        if( !mainbuf_.Load( *is_ ) ) return CResult( 0 );
    }

    return mainbuf_.ReadNext();
}

END_NS( srprism )
END_STD_SCOPES

