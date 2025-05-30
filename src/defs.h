#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RUNNING COMPILATION PROPERTIES :: COMMENT OR UNCOMMENT 

#define PROGRESS // DISPLAY % OF PROGRESS
#define PROFILE  // DISPLAY A PROFILE WITH EXISTS OR NOT

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef uint64_t U64;
typedef uint32_t U32;
typedef uint16_t U16;
typedef uint8_t  U8; 
typedef int64_t  I64;
typedef int32_t  I32;
typedef int16_t  I16;
typedef int8_t   I8;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define PNAME                  "AltaiR"
#define VERSION                1
#define RELEASE                2

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define MAX_STR                2048
#define BUFFER_SIZE            65535    
#define PROGRESS_MIN           200
#define DEF_HELP               0
#define DEF_FORCE              0
#define DEF_VERBOSE            0
#define DEF_VV                 0
#define DEF_IR                 1
#define DEF_AA                 1
#define DEF_OUT                0
#define DEF_PROF_G             1
#define DEF_PROF_R             1
#define DEF_SPLIT              1
#define DEF_PLOTS              0
#define DEF_MIN_CTX            11
#define DEF_MAX_CTX            15
#define DEF_THREADS            1
#define BGUARD                 32
#define ALPHABET_SIZE          4

#define MIN_NPARAM_FOR_PROGS   2

#define TARGET                 0
#define REFERENCE              1

#define DEF_FC_HELP            0
#define DEF_FC_VERBOSE         0
#define DEF_FC_COMPLETE        0
#define DEF_FC_MINIMUM         1
#define DEF_FC_MAXIMUM         99999999
#define DEF_FC_ALPHABET        "ACGT"
#define DEF_FC_SPECIAL         0

#define DEF_AF_HELP            0
#define DEF_AF_FL              0
#define DEF_AF_VERBOSE         0
#define DEF_AF_ALPHABET        "ACGT"

#define DEF_MA_HELP            0
#define DEF_MA_WINDOW          47
#define DEF_MA_COLUMN          1
#define DEF_MA_VERBOSE         0
#define DEF_MA_IGNORE          0
#define DEF_MA_POSITION        0

#define DEF_LR_HELP            0
#define DEF_LR_FORCE           0
#define DEF_LR_VERBOSE         0
#define DEF_LR_THRESHOLD       0
#define DEF_LR_WEIGHT          0.03
#define DEF_LR_IGNORE          10
#define DEF_LR_REGION          500
#define DEF_LR_ALPHABET        "ACGT"
#define DEF_LR_DNA             0
#define DEF_LR_MIN_LEVEL       1
#define DEF_LR_MAX_LEVEL       15
#define DEF_LR_LEVEL           11

#define DEF_NC_HELP            0
#define DEF_NC_FORCE           0
#define DEF_NC_VERBOSE         0
#define DEF_NC_THREADS         4
#define DEF_NC_ALPHABET        "ACGT"
#define DEF_NC_DNA             1
#define DEF_NC_MIN_LEVEL       1
#define DEF_NC_MAX_LEVEL       15
#define DEF_NC_LEVEL           11
#define DEF_NC_MIN_THREADS     1
#define DEF_NC_MAX_THREADS     999999
#define HEADERS_PREFIX_SIZE    50

#define DEF_NCD_HELP           0
#define DEF_NCD_FORCE          0
#define DEF_NCD_VERBOSE        0
#define DEF_NCD_THREADS        4
#define DEF_NCD_ALPHABET       "ACGT"
#define DEF_NCD_DNA            1
#define DEF_NCD_MIN_LEVEL      1
#define DEF_NCD_MAX_LEVEL      15
#define DEF_NCD_LEVEL          11
#define DEF_NCD_MIN_THREADS    1
#define DEF_NCD_MAX_THREADS    999999

#define DEF_SI_HELP            0
#define DEF_SI_FORCE           0
#define DEF_SI_VERBOSE         0
#define DEF_SI_ALPHABET        "ACGT"
#define DEF_SI_DNA             1

#define DEF_RW_HELP            0
#define DEF_RW_FORCE           0
#define DEF_RW_VERBOSE         0
#define DEF_RW_VV              0
#define DEF_RW_AA              1
#define DEF_RW_IR              1
#define DEF_RW_THREADS         0
#define DEF_RW_GC_P            0
#define DEF_RW_RW_P            0
#define DEF_RW_OUT             0
#define DEF_RW_PLOTS           0
#define DEF_RW_MIN             11
#define DEF_RW_MAX             14
#define DEF_RW_MIN_CTX         1
#define DEF_RW_MAX_CTX         20

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

