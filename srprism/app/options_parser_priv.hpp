#include <cassert>
#include <sstream>
#include <stdexcept>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
template< typename t_target >
class COptionsParser::CDefaultParser
{
    typedef COptionsParser::TValue TValue;

    public:

        void operator()( const TValue & val, t_target & target ) const
        {
            std::istringstream is( static_cast< const std::string & >( val ) );
            is >> target;
        }
};

template<>
class COptionsParser::CDefaultParser< std::string >
{
    typedef COptionsParser::TValue TValue;

    public:

        void operator()( const TValue & val, std::string & target ) const
        { target = static_cast< const std::string & >( val ); }
};

template<>
class COptionsParser::CDefaultParser< bool >
{
    typedef COptionsParser::TValue TValue;

    public:

        void operator()( const TValue & val, bool & target ) const
        {
            if( val == "true" ) target = true;
            else if( val == "false" ) target = false;
            else {
                throw std::domain_error( 
                        "value must be \"true\" or \"false\"" );
            }
        }
};

template<>
class COptionsParser::CDefaultParser< Uint1 >
{
    typedef COptionsParser::TValue TValue;

    public:

        void operator()( const TValue & val, Uint1 & target ) const
        {
            Uint4 tmp;
            std::istringstream is( static_cast< const std::string & >( val ) );
            is >> tmp;
            target = (Uint1)tmp;
        }
};

//------------------------------------------------------------------------------
namespace {
    enum EArgType { EVALUE = 0, EKEY, ESHORTKEY, ENONE };

    inline EArgType ArgType( const COptionsParser::TValue & arg )
    {
        if( arg.size() > 0 && arg[0] == '-' ) {
            if( arg.size() > 1 && arg[1] == '-' ) return EKEY;
            else if( arg.size() <= 1 || isalpha( arg[1] ) ) return ESHORTKEY;
        }
        
        return EVALUE;
    }

    inline const COptionsParser::TKey GetKey( 
            const COptionsParser::TValue & arg )
    { return arg.substr( 2 ); }

    inline const COptionsParser::TName GetShortKey(
            const COptionsParser::TValue & arg )
    { return arg.substr( 1 ); }
}

#define DCASENONE      state = ESTATE_END; break;
#define DCASEPOSVALUE  positionals_.push_back( *begin++ ); \
                       state = ESTATE_POS; \
                       break;
#define DCASESHORTKEY {\
                        TName short_key = GetShortKey( *begin++ ); \
                        TShortNameMap::iterator i_key = \
                            short_name_map_.find( short_key ); \
\
                        if( i_key == short_name_map_.end() ) { \
                            if( raise ) { \
                                M_THROW( CException, PARSE, \
                                        "unknown key " << short_key ); \
                            }\
                            else key = short_key;\
                        } \
                        else key = i_key->second; \
\
                        state = ProcessKey( key, raise ) ? \
                            ESTATE_VAL : ESTATE_KEY; \
                        break; \
                      }
#define DCASEKEY      key = GetKey( *begin++ ); \
                      state = ProcessKey( key, raise ) ? \
                        ESTATE_VAL : ESTATE_KEY; \
                      break;
#define DFAIL         default: SRPRISM_ASSERT( false );

template< typename t_fw_iter >
void COptionsParser::Parse( t_fw_iter begin, t_fw_iter end, bool raise )
{
    positionals_.clear();
    args_.clear();

    enum { ESTATE_START, 
           ESTATE_POS, 
           ESTATE_KEY, 
           ESTATE_VAL,
           ESTATE_END } state = ESTATE_START;
    TKey  key;

    do {
        EArgType arg_type = (begin == end) ? ENONE : ArgType( *begin );

        switch( state ) {
            case ESTATE_START: case ESTATE_POS:
                switch( arg_type ) {
                    case ENONE:     DCASENONE
                    case EVALUE:    DCASEPOSVALUE
                    case ESHORTKEY: DCASESHORTKEY
                    case EKEY:      DCASEKEY
                    DFAIL
                }

                break;

            case ESTATE_KEY:
                if( raise ) {
                    switch( arg_type ) {
                        case ENONE: case EKEY: case ESHORTKEY:

                            if( raise ) {
                                M_THROW( CException, PARSE, 
                                        "value expected after " << key );
                            }
                            else { state = ESTATE_END; break; }

                        case EVALUE: 
                            args_[key] = *begin++; state = ESTATE_VAL; break;
                        DFAIL
                    }
                }
                else {
                    switch( arg_type ) {
                        case ENONE:     state = ESTATE_END; break;
                        case EKEY:      DCASEKEY
                        case ESHORTKEY: DCASESHORTKEY
                        case EVALUE:
                            if( begin == end ) { state = ESTATE_END; break; }
                            else { 
                                args_[key] = *begin++; 
                                state = ESTATE_VAL; break; 
                            }
                        DFAIL
                    }
                }

                break;

            case ESTATE_VAL: 
                switch( arg_type ) {
                    case EVALUE: 

                        if( raise ) {
                            M_THROW( CException, PARSE,
                                    "expected option, got " << *begin );
                        }
                        else { state = ESTATE_END; break; }

                    case ENONE:     DCASENONE
                    case ESHORTKEY: DCASESHORTKEY
                    case EKEY:      DCASEKEY
                    DFAIL
                }

                break;

            DFAIL
        }
    } while( state != ESTATE_END );

    // check the number of positionals
    Uint4 minpos = optional_positionals_start_,
          maxpos = minpos + max_opt_positionals_;

    if( positionals_.size() < minpos ) {
        M_THROW( CException, MISSING, 
                 "too few positional arguments: expected at least " << minpos );
    }

    if( positionals_.size() > maxpos ) {
        M_THROW( CException, PARSE, 
                 "too many positional arguments: expected at most " << maxpos );
    }

    // assign default values to implicit positional parameters
    for( size_t i( positionals_.size() ); i < maxpos; ++i ) {
        positionals_.push_back( "" );
    }

    // assign defaults for and validate required parameters
    {
        typedef TKeys::const_iterator TKeyIter;
        typedef TParamMap::const_iterator TIter;

        for( TKeyIter i_key = required_keys_.begin(); 
                i_key != required_keys_.end(); ++i_key ) {
            if( args_.find( *i_key ) == args_.end() ) {
                TIter i_deflt = default_map_.find( *i_key );

                if( i_deflt != default_map_.end() ) {
                    args_[*i_key] = i_deflt->second;
                }
                else if( raise ){
                    M_THROW( CException, MISSING,
                             "missing required option " << *i_key );
                }
            }
        }
    }
}

#undef DFAIL
#undef DCASEKEY
#undef DCASESHORTKEY
#undef DCASEPOSVALUE
#undef DCASENONE

//------------------------------------------------------------------------------
template< typename t_target, typename t_target_parser >
void COptionsParser::Bind( 
        TParamPos pos, t_target & lval, const t_target_parser & target_parser )
{
    SRPRISM_ASSERT( pos < positionals_.size() );

    try {
        target_parser( positionals_[pos], lval );
    }
    catch( std::domain_error & e ) {
        M_THROW( CException, PARSE,
                 "error parsing value for positional argument " << pos <<
                 "; " << e.what() );
    }
}

//------------------------------------------------------------------------------
template< typename t_target, typename t_target_parser >
void COptionsParser::Bind( 
        const TKey & key, t_target & lval, 
        const t_target_parser & target_parser )
{
    SRPRISM_ASSERT( args_.find( key ) != args_.end() );

    try {
        target_parser( args_[key], lval );
    }
    catch( std::domain_error & e ) {
        M_THROW( CException, PARSE,
                 "error parsing value for option " << key <<
                 "; " << e.what() );
    }
}

END_NS( common )
END_STD_SCOPES

