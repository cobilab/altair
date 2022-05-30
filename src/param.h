#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <limits.h>

#include "defs.h"
#include "alphabet.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define NC_MAX_MEM                9999999
#define NC_MIN_MEM                0
#define NC_MAX_LEVEL              24
#define NC_MIN_LEVEL              1
#define NC_MAX_CTX                31
#define NC_MIN_CTX                1
#define NC_MAX_DEN                1000000
#define NC_MIN_DEN                1

#define NCD_MAX_MEM               9999999
#define NCD_MIN_MEM               0
#define NCD_MAX_LEVEL             24
#define NCD_MIN_LEVEL             1
#define NCD_MAX_CTX               31
#define NCD_MIN_CTX               1
#define NCD_MAX_DEN               1000000
#define NCD_MIN_DEN               1

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  char      **names;
  uint32_t  nFiles;
  }
SFILES;

typedef struct
  {
  U32       ctx;
  U32       den;
  U32       ir;
  U32       eIr;
  U32       memory;
  U32       hashSize;
  double    gamma;
  U32       edits;
  U32       eDen;
  double    eGamma;
  U8        type;
  U32       bet;
  double    subs_rate;
  double    adds_rate;
  double    dels_rate;
  U32       seed;
  U32       size;
  U32       init;
  U32       end;
  char      *name;
  }
MODEL_PAR;

typedef struct
  {
  U32       help;
  U32       verbose;
  U32       vv;
  U64       iBase;
  U32       force;
  U32       threads;
  U32       aa;
  U32       ir;
  U32       gc_prof;
  U32       raw_prof;
  U32       out;
  U32       plots;
  U32       min;
  U32       max;
  U32       nKmers;
  U32       split;
  U32       nSym;
  U64       *size;
  char      *parasite;
  char      *host;
  char      *prefix;
  }
RW_PARAMETERS;

typedef struct
  {
  ALPHABET  *A;
  char      *alphabet;
  U32       nSym;
  U8        help;
  U8        verbose;
  U8        force;
  U8        level;
  U8        dna;
  U32       threads;
  MODEL_PAR *model;
  U64       size;
  U32       nModels;
  U32       memory;
  char      *reference;
  char      *output;
  char      *filename;
  SFILES    *ref;
  }
NCD_PARAMETERS;

typedef struct
  {
  ALPHABET  *A;
  char      *alphabet;
  U32       nSym;
  U8        help;
  U8        verbose;
  U8        force;
  U8        level;
  U8        dna;
  U32       threads;
  MODEL_PAR *model;
  U64       size;
  U32       nModels;
  U32       memory;
  char      *reference;
  char      *output;
  char      *filename;
  SFILES    *ref;
  }
NC_PARAMETERS;

typedef struct
  {
  U32       help;
  U32       verbose;
  U32       ignore;
  U32       show_pos;
  U32       column;
  U32       window_size;
  }
MA_PARAMETERS;

typedef struct
  {
  U32       help;
  U32       verbose;
  U32       f_line;
  char      *alphabet;
  }
AF_PARAMETERS;

typedef struct
  {
  U32       help;
  U32       verbose;
  U32       complete;
  U32       minimum;
  U32       maximum;
  char      *alphabet;
  double    cg_min;
  double    cg_max;
  char      **patterns;
  U32       nPatterns;
  char      **ignore;
  U32       nIgnore;
  }
FC_PARAMETERS;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t    ArgsNum            (uint32_t , char *[], uint32_t, char *, char *,
                               uint32_t, uint32_t);
double      ArgsDouble         (double, char *[], uint32_t, char *, char *);
uint8_t     ArgsState          (uint8_t  , char *[], uint32_t, char *, char *);
char        *ArgsString        (char    *, char *[], uint32_t, char *, char *);
char        *ArgsFiles         (char *[], uint32_t, char *);
MODEL_PAR   ArgsUniqModelNC    (char *, uint8_t);
MODEL_PAR   ArgsUniqModelNCD   (char *, uint8_t);
void        PrintParametersMA  (MA_PARAMETERS *);
void        PrintParametersAF  (AF_PARAMETERS *);
void        PrintParametersFC  (FC_PARAMETERS *);
void        PrintParametersNC  (NC_PARAMETERS *);
void        PrintParametersNCD (NCD_PARAMETERS *);
void        PrintParametersRW  (RW_PARAMETERS *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
