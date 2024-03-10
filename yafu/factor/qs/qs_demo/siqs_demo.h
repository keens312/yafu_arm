/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

//definitions
#ifndef _SIQD_DEMO_H
#define _SIQD_DEMO_H

#define VERSION_STRING "1.1"

//default maximum size in chars for a str_t
#define GSTR_MAXSIZE 1024

//support libraries
#include <stdint.h>
#include <gmp.h>

#ifdef _MSC_VER
#include <Windows.h>
#endif

// a structure to hold a bunch of configuration info
// for yafu, instead of declaring a bunch of globals.
typedef struct
{
    // isprime
    uint32_t NUM_WITNESSES;

    // output behavior - used everywhere
    int VFLAG, LOGFLAG;

    // threading
    int THREADS;

    // input options
    int USEBATCHFILE;
    int USERSEED;
    int CMD_LINE_REPEAT;
    char batchfilename[1024];
    char sessionname[1024];
    char scriptname[1024];

    // machine info
    double MEAS_CPU_FREQUENCY;
    int VERBOSE_PROC_INFO;
    char CPU_ID_STR[80];
    uint32_t L1CACHE, L2CACHE;
    int CLSIZE;
    char HAS_SSE41;
    char HAS_AVX;
    char HAS_AVX2;
#if defined(WIN32)
    char sysname[MAX_COMPUTERNAME_LENGTH + 1];
    int sysname_sz;
#else
    char sysname[256];
    int sysname_sz;
#endif  

} siqs_obj_t;

extern void siqs_init(siqs_obj_t* obj);

#endif //ifndef HEAD_DEF
