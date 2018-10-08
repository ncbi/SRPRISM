#ifndef __AM_COMMON_JAILDEF_H__
#define __AM_COMMON_JAILDEF_H__

#include <config.h>

#if HAVE_INTTYPES_H == 1
#   include <inttypes.h>
#endif

//------------------------------------------------------------------------------
// namespace setup
//

#ifdef ROOT_NS
#   define START_ROOT_NS namespace ROOT_NS {
#   define END_ROOT_NS }
#   define USE_ROOT_NS using namespace ROOT_NS;
#else
#   define START_ROOT_NS
#   define END_ROOT_NS
#   define USE_ROOT_NS
#endif

#ifdef JAIL_NS
#   define START_JAIL_NS namespace JAIL_NS {
#   define END_JAIL_NS }
#   define USE_JAIL_NS using namespace JAIL_NS;
#else
#   define START_JAIL_NS
#   define END_JAIL_NS
#   define USE_JAIL_NS
#endif

#ifdef ROOT_NS
#   ifdef JAIL_NS
#       define STD_SCOPES ROOT_NS::JAIL_NS
#   else
#       define STD_SCOPES ROOT_NS
#   endif
#else
#   ifdef JAIL_NS
#       define STD_SCOPES JAIL_NS
#   else
#       define STD_SCOPES
#   endif
#endif

#define START_NS(ns) namespace ns {
#define END_NS(ns) }
#define USE_NS(ns) using namespace ns;

START_ROOT_NS
START_JAIL_NS
START_NS( common )

//------------------------------------------------------------------------------
// fixed width integer types
//

#if HAVE_INTTYPES_H == 1

typedef uint8_t  Uint1;
typedef uint16_t Uint2;
typedef uint32_t Uint4;
typedef uint64_t Uint8;

typedef int8_t  Sint1;
typedef int16_t Sint2;
typedef int32_t Sint4;
typedef int64_t Sint8;

#elif defined(WIN32)

typedef unsigned __int8  Uint1;
typedef unsigned __int16 Uint2;
typedef unsigned __int32 Uint4;
typedef unsigned __int64 Uint8;

typedef signed __int8  Sint1;
typedef signed __int16 Sint2;
typedef signed __int32 Sint4;
typedef signed __int64 Sint8;

#endif

END_NS( common )
END_JAIL_NS
END_ROOT_NS

#endif

