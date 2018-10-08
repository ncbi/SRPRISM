/*  $Id: mkidx.cpp 426095 2014-02-05 20:58:39Z morgulis $
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
 * File Description: index creation code
 *
 */

#include <ncbi_pch.hpp>

#include <memory>

#include "../common/util.hpp"
#include "../common/textfile.hpp"
#include "../seq/seqinput_factory.hpp"
#include "../seq/seqinput.hpp"
#include "memmgr.hpp"
#include "seqstore.hpp"
#include "mkidx.hpp"
#include "seqstore_factory.hpp"
#include "mkidx_pass.hpp"
#include "nmer_iterator.hpp"

START_STD_SCOPES
START_NS( srprism )
USE_NS( common )
USE_NS( seq )

//------------------------------------------------------------------------------
namespace {
    const char * INFMT_FASTA = "fasta";
    const char * INFMT_FASTQ = "fastq";

    const char * OUTFMT_STANDARD = "standard";

    typedef TSeqSize TNPos;

    struct SPosTableEntry
    {
        CSeqStore::TPos pos;
        TPrefix prefix;
        TStrand strand;
    };
}

//------------------------------------------------------------------------------
void CMkIdx::SaveIdxHeader( CWriteBinFile & idx_file )
{
    char e( (char)Endianness() );
    idx_file.Write( &e, 1 );
    e = (char)(PREFIX_LEN - MAP_PREFIX_LEN);
    idx_file.Write( &e, 1 );
    e = 0;
    idx_file.Write( &e, 1 );
    idx_file.Write( &e, 1 );
}

//------------------------------------------------------------------------------
CMkIdx::CMkIdx( const SOptions & options )
    : alt_loc_spec_name_( options.alt_loc_spec ),
      output_( options.output ),
      infmt_(  options.infmt ),
      outfmt_( options.outfmt ),
      input_c_( options.input_compression ),
      max_mem_( MEGABYTE*options.max_mem ),
      ss_seg_len_( options.ss_seg_len ),
      al_extend_( options.al_extend )
{
    if( !options.input.empty() ) {
        std::string::size_type pos( 0 ), pos1;

        do {
            pos1 = options.input.find_first_of( ",", pos );
            input_.push_back( options.input.substr( pos, pos1 - pos ) );
            pos = pos1 + 1;
        } while( pos1 != std::string::npos );
    }
    else if( !options.input_list.empty() ) {
        std::auto_ptr< CReadTextFile > lst( 
                CReadTextFile::MakeReadTextFile( options.input_list ) );

        while( !lst->Eof() ) {
            std::string line( lst->GetLine() );
            if( line.empty() || line[0] == '#' ) continue;
            input_.push_back( line );
        }
    }
    else input_.push_back( "" );

    Validate();
}

//------------------------------------------------------------------------------
void CMkIdx::Validate( void )
{
    { // validation of input format
        static const char * FMT_NAMES[] = { INFMT_FASTA, INFMT_FASTQ, 0 };
        bool found( false );

        for( const char **fmt_name = FMT_NAMES; *fmt_name != 0; ++fmt_name ) {
            if( infmt_ == *fmt_name ) { found = true; break; }
        }

        if( !found ) {
            std::string valid_formats = "valid formats:";

            for( const char **fmt_name = FMT_NAMES; 
                    *fmt_name != 0; ++fmt_name ) {
                valid_formats += " ";
                valid_formats += *fmt_name;
            }

            M_THROW( CException, VALIDATE,
                     "input format " << infmt_ << " is not recognized." <<
                     "(" << valid_formats << ")" );
        }
    }

    { // validation of output format
        static const char * FMT_NAMES[] = { OUTFMT_STANDARD, 0 };
        bool found( false );

        for( const char **fmt_name = FMT_NAMES; *fmt_name != 0; ++fmt_name ) {
            if( outfmt_ == *fmt_name ) { found = true; break; }
        }

        if( !found ) {
            std::string valid_formats = "valid formats:";

            for( const char **fmt_name = FMT_NAMES; 
                    *fmt_name != 0; ++fmt_name ) {
                valid_formats += " ";
                valid_formats += *fmt_name;
            }

            M_THROW( CException, VALIDATE,
                     "output format " << outfmt_ << " is not recognized." <<
                     "(" << valid_formats << ")" );
        }
    }

    { // validation of output base name
        if( output_.empty() ) {
            M_THROW( CException, VALIDATE,
                     "output base name should not be empty" );
        }
    }

    { // validation of memory limit
        if( max_mem_ == 0 ) {
            M_THROW( CException, VALIDATE,
                     "max memory must be positive" );
        }
    }
}

//------------------------------------------------------------------------------
void CMkIdx::MkSeqStore( void )
{
    M_TRACE( CTracer::INFO_LVL, "creating sequence store" );
    CSeqStoreFactory seqstore( 
            max_mem_, output_, alt_loc_spec_name_, ss_seg_len_, al_extend_ );

    for( std::vector< std::string >::const_iterator ii( input_.begin() );
            ii != input_.end(); ++ii ) {
        std::auto_ptr< CSeqInput > seq_in( CSeqInputFactory::MakeSeqInput( 
                    infmt_, *ii, 1, input_c_ ) );

        while( !seq_in->Done() ) {
            if( !seq_in->Next() ) break;
            seqstore.Append( seq_in->Id(), seq_in->Data( 0 ) );
        }
    }

    seqstore.Save();
}

//------------------------------------------------------------------------------
void CMkIdx::Run( void )
{
    MkSeqStore();

    CMemoryManager mem_mgr( max_mem_ );
    CSeqStore seq_store( output_, mem_mgr );
    seq_store.Load();
    
    static const TSeqSize HASH_KEY_SIZE = CMkIdxPass::HASH_KEY_SIZE;
    static const size_t HASH_KEY_BITS = 
        HASH_KEY_SIZE*SCodingTraits< SEQDATA_CODING >::LETTER_BITS;
    static const size_t SHIFT = BYTEBITS*sizeof( TPrefix ) - HASH_KEY_BITS;
    static const size_t NUM_HASH_KEYS = 
        SBitFieldTraits< size_t, HASH_KEY_BITS >::MAX + 1;

    size_t * counts_table( 
            (size_t *)mem_mgr.Allocate( NUM_HASH_KEYS*sizeof( size_t ) ) );
    std::fill( counts_table, counts_table + NUM_HASH_KEYS, (size_t)0 );

    {
        M_TRACE( CTracer::INFO_LVL, "generating n-mer counts" );
        
        for( size_t seq_idx = 0; seq_idx < seq_store.NSeq(); ++seq_idx ) {
            for( CNMerIterator nmer_iter( seq_store, seq_idx ); 
                    !nmer_iter.End(); nmer_iter.Next() ) {
                ++counts_table[nmer_iter.Prefix()>>SHIFT];
            }
        }
    }

    CWriteBinFile idx_file( output_ + IDX_PROPER_SFX );
    SaveIdxHeader( idx_file );
    CWriteBinFile map_file( output_ + IDX_MAP_SFX );
    CWriteBinFile rmap_file( output_ + IDX_REPMAP_SFX );
    size_t hash_key_start( 0 );
    size_t free_space_size( mem_mgr.GetFreeSpace() );
    void * free_space( mem_mgr.Allocate( free_space_size ) );
    size_t pass_no( 1 );
    M_TRACE( CTracer::INFO_LVL, "generating index" );

    while( hash_key_start < NUM_HASH_KEYS ) {
        M_TRACE( CTracer::INFO_LVL, "PASS " << pass_no );
        size_t t( hash_key_start );
        CMkIdxPass pass( 
                counts_table, seq_store, free_space, free_space_size,
                hash_key_start, map_file, rmap_file, idx_file );
        if( hash_key_start == t ) M_THROW( CException, MEMORY, "" );
        pass.Run();
        ++pass_no;
    }

    mem_mgr.Free( free_space );
    mem_mgr.Free( (void *)counts_table );
    seq_store.Unload();
}

END_NS( srprism )
END_STD_SCOPES

