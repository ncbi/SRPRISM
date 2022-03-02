#ifndef __AM_COMMON_DEF_H__
#define __AM_COMMON_DEF_H__

#include <cstdlib>
#include <cassert>

/*
#ifdef WIN32
#	ifndef NCBI_CPP_TK
#		define NCBI_CPP_TK 1
#	endif
#endif
*/

/*
#ifndef NCBI_CPP_TK

#include <common/jaildef.h>

#else

#include <../src/internal/align_toolbox/srprism/lib/common/jaildef.h>

#endif
*/
#include "../common/jaildef.h"

//------------------------------------------------------------------------------
// namespace setup
//

#define USE_STD_SCOPES   USE_ROOT_NS USE_JAIL_NS
#define START_STD_SCOPES START_ROOT_NS START_JAIL_NS
#define END_STD_SCOPES   END_JAIL_NS END_ROOT_NS

START_STD_SCOPES
START_NS( common )

//------------------------------------------------------------------------------
// system dependent info
//
const char FPATH_SEP = '/';

//------------------------------------------------------------------------------
// some useful constants
//
const size_t BYTEBITS = 8;
const size_t KILOBYTE = 1024ULL;
const size_t MEGABYTE = 1024*KILOBYTE;

//------------------------------------------------------------------------------
// toolkit dependent definitions
//
#ifdef NCBI_CPP_TK
#   define SRPRISM_ASSERT _ASSERT
#else
#   define SRPRISM_ASSERT assert
#endif

END_NS( common )
END_STD_SCOPES

#endif
