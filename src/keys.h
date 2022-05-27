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
#define K9  9

typedef struct
  {
  char *key;
  int val;
  }
K_STRUCT;

static K_STRUCT LT_KEYS[] = 
  { 
    { "help"        , K1  },  // HELP MENU
    { "filter"      , K2  },  // FILTER CHARACTERISTICS
    { "frequency"   , K3  },  // ALPHABET FREQUENCIES
    { "average"     , K4  },  // MOVING AVERAGE
    { "redundancy"  , K5  },  // REDUNDANCY PROFILES
    { "nc"          , K6  },  // NORMALIZED COMPRESSION
    { "ncd"         , K7  },  // NORMALIZED COMPRESSION DISTANCE
    { "simulate"    , K8  },  // SIMULATION OF FASTA CONTENT
    { "raw"         , K9  }   // RELATIVE ABSENT WORDS
  };

#define NKEYS (sizeof(LT_KEYS)/sizeof(K_STRUCT))

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int KeyString (char *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

