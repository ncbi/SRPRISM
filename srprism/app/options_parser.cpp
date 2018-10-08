#include <common/def.h>

#include <string>
#include <sstream>
#include <algorithm>

#include <common/exception.hpp>
#include <common/text_formatter.hpp>
#include <options_parser.hpp>

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
const COptionsParser::TValue COptionsParser::TRUE_VAL =  "true", 
                             COptionsParser::FALSE_VAL = "false";

//------------------------------------------------------------------------------
COptionsParser::COptionsParser( 
        const std::string program_name,
        const std::string & program_description,
        const std::string & version )
    : optional_positionals_start_( 0 ), 
      prog_descr_( program_description + "\nversion " + version ),
      max_opt_positionals_( 0 )
{
    synopsis_ += std::string( "USAGE:\n\n\t" ) + program_name;
}

//------------------------------------------------------------------------------
const std::string COptionsParser::Usage() const
{
    std::ostringstream os;
    os << "\n" 
       << CTextFormatter( 
               FMT_PROG_DESCR_START, 
               FMT_LINE_RIGHT - FMT_PROG_DESCR_START )( prog_descr_ )
        << "\n\n";
    os << synopsis_ << "\n\n";
    
    {
        typedef TGroups::const_iterator TIter;

        for( TIter i = groups_.begin(); i != groups_.end(); ++i ) {
            os << i->header << "\n\n" << i->data;
        }
    }

    os << std::endl;
    return os.str();
}

//------------------------------------------------------------------------------
void COptionsParser::NewGroup( const std::string & title )
{
    SHelpGroup newgroup = { title };
    groups_.push_back( newgroup );
}

//------------------------------------------------------------------------------
void COptionsParser::AddPositionalDescription(
        const std::string & param_label, 
        const std::string & description )
{
    if( groups_.empty() ) NewGroup( "UNKNOWN" );
    std::string & help = groups_.rbegin()->data;

    if( required_keys_.empty() && optional_keys_.empty() && flags_.empty() ) {
        synopsis_ += " " + param_label;
        help.append( 
                CTextFormatter( 
                    FMT_PARAM_START, 
                    FMT_LINE_RIGHT - FMT_PARAM_START 
                )( param_label + " [required]\n" ) );
        help.append(
                CTextFormatter(
                    FMT_PARAM_DESCR_START,
                    FMT_LINE_RIGHT - FMT_PARAM_DESCR_START
                )( description ) );
        ++optional_positionals_start_;
    }
    else {
        M_THROW( CException, ORDER,
                 "positional parameter appears after an option" );
    }
}

//------------------------------------------------------------------------------
void COptionsParser::AddParamPriv( 
        const TKey & key, const TName & short_name, TKeys & keys )
{
    if( std::find( keys.begin(), keys.end(), key ) != keys.end() ) {
        M_THROW( CException, KEY, "duplicate key " << key );
    }

    if( !short_name.empty() && 
        short_name_map_.find( short_name ) != short_name_map_.end() ) {
        M_THROW( CException, KEY, "duplicate short key " << short_name );
    }

    keys.push_back( key );
    if( !short_name.empty() ) short_name_map_[short_name] = key;
}

//------------------------------------------------------------------------------
void COptionsParser::AddFlag( 
        const TKey & key, const TName & short_name, 
        const std::string & description )
{
    AddParamPriv( key, short_name, required_keys_ );
    flags_.push_back( key );
    default_map_[key] = FALSE_VAL;
    synopsis_ += " [--"  + key + "]";
    std::string suffix = " [flag]\n";

    if( !short_name.empty() ) {
        suffix = std::string( "|-" ) + short_name + suffix;
    }

    if( groups_.empty() ) NewGroup( "UNKNOWN" );
    std::string & help = groups_.rbegin()->data;

    help.append( 
            CTextFormatter( 
                FMT_PARAM_START, 
                FMT_LINE_RIGHT - FMT_PARAM_START 
            )( std::string( "--" ) + key + suffix ) );
    help.append( 
            CTextFormatter( 
                FMT_PARAM_DESCR_START,
                FMT_LINE_RIGHT - FMT_PARAM_DESCR_START
            )( description ) );
}

//------------------------------------------------------------------------------
void COptionsParser::AddParam( 
        const TKey & key, const TName & short_name, 
        const std::string & description, const std::string & param_label )
{
    AddParamPriv( key, short_name, required_keys_ );
    synopsis_ += " [--" + key + " <" + param_label + ">]";
    std::string suffix = std::string( " <" ) + param_label + "> [required]\n";

    if( !short_name.empty() ) {
        suffix = std::string( "|-" ) + short_name + suffix;
    }

    if( groups_.empty() ) NewGroup( "UNKNOWN" );
    std::string & help = groups_.rbegin()->data;

    help.append( 
            CTextFormatter( 
                FMT_PARAM_START, 
                FMT_LINE_RIGHT - FMT_PARAM_START 
            )( std::string( "--" ) + key + suffix ) );
    help.append( 
            CTextFormatter( 
                FMT_PARAM_DESCR_START,
                FMT_LINE_RIGHT - FMT_PARAM_DESCR_START
            )( description ) );
}

//------------------------------------------------------------------------------
void COptionsParser::AddDefaultParam( 
        const TKey & key, const TName & short_name, const TValue & default_val,
        const std::string & description, const std::string & param_label )
{
    AddParamPriv( key, short_name, required_keys_ );
    default_map_[key] = default_val;
    synopsis_ += " [--" + key + " <" + param_label + ">]";
    std::string suffix = std::string( " <" ) + param_label
                       + "> [default: " + default_val + "]\n";

    if( !short_name.empty() ) {
        suffix = std::string( "|-" ) + short_name + suffix;
    }

    if( groups_.empty() ) NewGroup( "UNKNOWN" );
    std::string & help = groups_.rbegin()->data;

    help.append( 
            CTextFormatter( 
                FMT_PARAM_START, 
                FMT_LINE_RIGHT - FMT_PARAM_START 
            )( std::string( "--" ) + key + suffix ) );
    help.append( 
            CTextFormatter( 
                FMT_PARAM_DESCR_START,
                FMT_LINE_RIGHT - FMT_PARAM_DESCR_START
            )( description ) );
}

//------------------------------------------------------------------------------
void COptionsParser::AddOptionalParam( 
        const TKey & key, const TName & short_name, 
        const std::string & description, const std::string & param_label )
{
    AddParamPriv( key, short_name, optional_keys_ );
    synopsis_ += " [--" + key + " <" + param_label + ">]";
    std::string suffix = std::string( " <" ) + param_label + "> [optional]\n";

    if( !short_name.empty() ) {
        suffix = std::string( "|-" ) + short_name + suffix;
    }

    if( groups_.empty() ) NewGroup( "UNKNOWN" );
    std::string & help = groups_.rbegin()->data;

    help.append( 
            CTextFormatter( 
                FMT_PARAM_START, 
                FMT_LINE_RIGHT - FMT_PARAM_START 
            )( std::string( "--" ) + key + suffix ) );
    help.append( 
            CTextFormatter( 
                FMT_PARAM_DESCR_START,
                FMT_LINE_RIGHT - FMT_PARAM_DESCR_START
            )( description ) );
}

//------------------------------------------------------------------------------
bool COptionsParser::ProcessKey( const TKey & key, bool raise )
{
    bool flag = std::find( flags_.begin(), flags_.end(), key ) != flags_.end();

    if( raise ) {
        bool found =
            flag ||
            std::find( required_keys_.begin(), required_keys_.end(), key ) !=
                required_keys_.end() ||
            std::find( optional_keys_.begin(), optional_keys_.end(), key ) !=
                optional_keys_.end();
        if( !found ) M_THROW( CException, KEY, "unknown key " << key );
    }

    if( flag ) args_[key] = TRUE_VAL;
    return flag;
}

END_NS( common )
END_STD_SCOPES

