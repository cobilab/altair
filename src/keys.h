#ifndef KEYS_H_INCLUDED
#define KEYS_H_INCLUDED

#include <stdio.h>
#include <stdint.h>

#include "defs.h"

#define K0  0
#define K1  1
#define K2  2
#define K3  3
#define K4  4
#define K5  5
#define K6  6
#define K7  7
#define K8  8

typedef struct
  {
  char *key;
  int val;
  }
K_STRUCT;

static K_STRUCT LT_KEYS[] = 
  { 
    { "help"        , K1  },  // HELP MENU
    { "average"     , K2  },  // MOVING AVERAGE
    { "filter"      , K3  },  // FILTER CHARACTERISTICS
    { "frequency"   , K4  },  // ALPHABET FREQUENCIES
    { "nc"          , K5  },  // NORMALIZED COMPRESSION
    { "ncd"         , K6  },  // NORMALIZED COMPRESSION DISTANCE
    { "raw"         , K7  }   // RELATIVE ABSENT WORDS
  };

#define NKEYS (sizeof(LT_KEYS)/sizeof(K_STRUCT))

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int KeyString (char *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

