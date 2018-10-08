#ifndef __AM_COMMON_OPTIONS_PARSER_HPP__
#define __AM_COMMON_OPTIONS_PARSER_HPP__

#include <common/def.h>

#include <string>
#include <vector>
#include <map>

#include <common/exception.hpp>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
class COptionsParser
{
    private:

        template< typename t_target > struct CDefaultParser;

    public:

        struct CException : public common::CException
        {
            typedef common::CException TBase;

            static const TErrorCode ORDER   = 0;
            static const TErrorCode PARSE   = 1;
            static const TErrorCode KEY     = 2;
            static const TErrorCode MISSING = 3;

            virtual const std::string ErrorMessage( TErrorCode code ) const
            {
                if( code == ORDER ) return "order error";
                else if( code == PARSE ) return "options parser error";
                else if( code == KEY ) return "key error";
                else if( code == MISSING ) return "argument missing";
                else return TBase::ErrorMessage( code );
            }

            M_EXCEPT_CTOR( CException )
        };

        typedef std::string TName;
        typedef std::string TValue;
        typedef TName TKey;
        typedef Uint2 TParamPos;

        static const TValue TRUE_VAL, FALSE_VAL;

        COptionsParser( 
                const std::string program_name,
                const std::string & program_description,
                const std::string & version = "" );

        template< typename t_fw_iter >
        void Parse( t_fw_iter begin, t_fw_iter end, bool raise = true );

        void AddParam( 
                const TKey & key, 
                const TName & short_name,
                const std::string & description,
                const std::string & param_label = "arg" );

        void AddOptionalParam( 
                const TKey & key, 
                const TName & short_name,
                const std::string & description,
                const std::string & param_label = "arg" );

        void AddFlag( 
                const TKey & key, 
                const TName & short_name,
                const std::string & description );

        void AddDefaultParam( 
                const TKey & key, 
                const TName & short_name,
                const TValue & default_val,
                const std::string & description,
                const std::string & param_label = "arg" );

        void AddPositionalDescription(
                const std::string & param_label,
                const std::string & description );

        void AddPositionalDescriptionWithDefault(
                const TValue & default_val,
                const std::string & param_label,
                const std::string & description );

        bool IsPresent( const TKey & key ) const
        { return args_.find( key ) != args_.end(); }

        template< typename t_target, typename t_target_parser >
        void Bind( 
                const TKey & key,
                t_target & lval, 
                const t_target_parser & target_parser );

        template< typename t_target, typename t_target_parser >
        void Bind( 
                TParamPos pos,
                t_target & lval, 
                const t_target_parser & target_parser );

        template< typename t_target > 
        void Bind( const TKey & key, t_target & lval )
        { Bind( key, lval, COptionsParser::CDefaultParser< t_target >() ); }

        template< typename t_target > 
        void Bind( TParamPos pos, t_target & lval )
        { Bind( pos, lval, COptionsParser::CDefaultParser< t_target >() ); }

        const std::string Usage() const;
        TParamPos NPositionals() const { return positionals_.size(); }

        void SetMaxOptPositionals( Uint4 val ) { 
            max_opt_positionals_ = val; 
        }

        void NewGroup( const std::string & title );

    private:

        typedef std::vector< TValue > TPositionals;
        typedef std::vector< TKey > TKeys;
        typedef std::map< TName, TKey > TShortNameMap;
        typedef std::map< TKey, TValue > TParamMap;

        static const Uint1 FMT_PROG_DESCR_START = 0;
        static const Uint1 FMT_PARAM_START = 8;
        static const Uint1 FMT_PARAM_DESCR_START = 8;
        static const Uint1 FMT_LINE_RIGHT = 70;

        struct SHelpGroup
        {
            std::string header;
            std::string data;
        };

        typedef std::vector< SHelpGroup > TGroups;

        bool ProcessKey( const TKey & key, bool raise );

        void AddParamPriv( 
                const TKey & key, 
                const TName & short_name,
                TKeys & keys );

        TPositionals positionals_;
        TParamPos optional_positionals_start_;
        TKeys required_keys_;
        TKeys optional_keys_;
        TKeys flags_;
        TShortNameMap short_name_map_;
        TParamMap default_map_;
        TParamMap args_;
        std::string synopsis_;
        const std::string prog_descr_;
        TGroups groups_;
        Uint4 max_opt_positionals_;
};

END_NS( common )
END_STD_SCOPES

#include <options_parser_priv.hpp>
#endif

