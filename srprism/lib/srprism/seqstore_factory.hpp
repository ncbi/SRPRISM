/*  $Id: seqstore_factory.hpp 637057 2021-09-05 23:00:51Z morgulis $
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
 * File Description: seqstore generation
 *
 */

#ifndef __SRPRISM_SEQSTORE_FACTORY_HPP__
#define __SRPRISM_SEQSTORE_FACTORY_HPP__

#include "../common/def.h"

#include <string>
#include <vector>

#include "../common/exception.hpp"
#include "../common/textfile.hpp"
#include "../common/binfile.hpp"
#include "../seq/seqdef.hpp"
#include "memmgr.hpp"
#include "seqstore_base.hpp"

START_STD_SCOPES
START_NS( srprism )

//------------------------------------------------------------------------------
class CSeqStoreFactory : public CSeqStoreBase
{
    private:

        typedef std::vector< std::string > TIdMap;
        typedef std::vector< seq::TLetter  > TAmbigData;
        typedef std::vector< TWord * > TSeqData;

        struct SAmbigRunData
        {
            TDBOrdId sid;
            seq::TSeqSize pos;
            seq::TSeqSize len;
            size_t offset;

            SAmbigRunData( TDBOrdId s, seq::TSeqSize p, size_t o ) 
                : sid( s ), pos( p ), len( 0 ), offset( o )
            {}

            friend bool operator<( 
                    const SAmbigRunData & l, const SAmbigRunData & r )
            { 
                if( l.sid == r.sid ) return l.pos < r.pos;
                return l.sid < r.sid;
            }
        };

        typedef std::vector< SAmbigRunData > TAmbigMap;

        struct SRevIdMapEntry
        {
            std::string id;
            TDBOrdId oid;

            SRevIdMapEntry( const std::string & arg_id, TDBOrdId arg_oid )
                : id( arg_id ), oid( arg_oid )
            {}

            friend bool operator<( 
                    const SRevIdMapEntry & l, const SRevIdMapEntry & r )
            { return l.id < r.id; }
        };

        typedef std::vector< SRevIdMapEntry > TRevIdMap;

        struct SSeqInfoEntry
        {
            TDBOrdId oid;
            TDBOrdId ref_oid;
            seq::TSeqSize alt_loc_start,
                          alt_loc_end,
                          alt_loc_real_end,
                          ref_loc_start,
                          ref_loc_end;
            bool fuzzy_left,
                 fuzzy_right,
                 flip;

            SSeqInfoEntry( TDBOrdId arg_oid, TSeqSize len )
                : oid( arg_oid ), ref_oid( oid ),
                  alt_loc_start( 0 ), alt_loc_end( len ),
                  alt_loc_real_end( len ),
                  ref_loc_start( 0 ), ref_loc_end( len ),
                  fuzzy_left( false ), fuzzy_right( false ), flip( false )
            {}

            SSeqInfoEntry( 
                    TDBOrdId arg_oid, TDBOrdId arg_ref_oid,
                    seq::TSeqSize arg_alt_loc_start,
                    seq::TSeqSize arg_alt_loc_end,
                    seq::TSeqSize arg_alt_loc_real_end,
                    seq::TSeqSize arg_ref_loc_start,
                    seq::TSeqSize arg_ref_loc_end,
                    bool arg_fuzzy_left, bool arg_fuzzy_right, bool arg_flip )
                : oid( arg_oid ), ref_oid( arg_ref_oid ),
                  alt_loc_start( arg_alt_loc_start ),
                  alt_loc_end( arg_alt_loc_end ),
                  alt_loc_real_end( arg_alt_loc_real_end ),
                  ref_loc_start( arg_ref_loc_start ),
                  ref_loc_end( arg_ref_loc_end ),
                  fuzzy_left( arg_fuzzy_left ),
                  fuzzy_right( arg_fuzzy_right ),
                  flip( arg_flip )
            {}

            friend bool operator<(
                    const SSeqInfoEntry & l, const SSeqInfoEntry & r )
            { return l.oid < r.oid; }
        };

        typedef std::vector< SSeqInfoEntry > TSeqInfo;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode MEMORY = 0;
            static const TErrorCode ALTLOC = 1;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == MEMORY ) return "memory limit exceeded";
                else if( code == ALTLOC ) {
                    return "alternative loci specification error";
                }
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        CSeqStoreFactory( 
                size_t max_mem, 
                const std::string & base_name,
                const std::string & alt_loc_spec_name,
                size_t segment_letters,
                size_t al_extend );

        template< typename data_t >
        void Append( const std::string & id, const data_t & seq_data );

        void Save( void );

    private:

        static std::string GetWord( 
                const std::string & line, 
                std::string::size_type & pos );

        void ComputeSeqLetters( void );
        void SaveHeader( void );
        void SetUpSeqInfo( void );
        void ParseAltLocLine( const std::string & line );
        void SaveSeqData( TDBOrdId oid );

        static void StoreSeqSeg( 
                const TWord * seg_data, TSeqSize seg_start, TSeqSize seg_len,
                TWord * dest, TSeqSize dest_start );

        void StoreAmbigInfo(
                TDBOrdId oid, TSeqSize start_pos, TSeqSize len, TSeqSize offset,
                std::vector< SAmbigRun > & amap, 
                std::vector< TLetter > & adata );

        void StoreAmbigInfoReverse(
                TDBOrdId oid, TSeqSize start_pos, TSeqSize len, TSeqSize offset,
                std::vector< SAmbigRun > & amap, 
                std::vector< TLetter > & adata );

        CMemoryManager mem_mgr_;
        std::string base_name_;
        std::string alt_loc_spec_name_;
        TIdMap id_map_;
        TRevIdMap rev_id_map_;
        TSeqInfo seq_info_;
        TAmbigMap ambig_map_;
        TAmbigData ambig_data_;
        TAmbigMask ambig_mask_;
        TSeqData seq_data_,
                 mask_data_;
        std::unique_ptr< common::CWriteBinFile > ss_outs_;
        std::unique_ptr< common::CWriteBinFile > mss_outs_;
        std::unique_ptr< common::CWriteTextFile > idmap_outs_;
        std::unique_ptr< common::CWriteBinFile > seqmap_outs_;
        std::unique_ptr< common::CWriteBinFile > ambig_map_outs_;
        std::unique_ptr< common::CWriteBinFile > ambig_data_outs_;
        /*
        std::auto_ptr< common::CWriteBinFile > ss_outs_;
        std::auto_ptr< common::CWriteBinFile > mss_outs_;
        std::auto_ptr< common::CWriteTextFile > idmap_outs_;
        std::auto_ptr< common::CWriteBinFile > seqmap_outs_;
        std::auto_ptr< common::CWriteBinFile > ambig_map_outs_;
        std::auto_ptr< common::CWriteBinFile > ambig_data_outs_;
        */
        TPos curr_pos_;
        size_t ambig_map_size_;
        size_t ambig_data_size_;
        size_t segment_letters_;
        size_t min_seq_len_;
        TSeqSize al_extend_;
};

END_NS( srprism )
END_STD_SCOPES

#include "seqstore_factory_priv.hpp"

#endif

