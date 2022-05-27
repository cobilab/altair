#ifndef FC_H_INCLUDED
#define FC_H_INCLUDED

#include <stdio.h>
#include <stdint.h>

#include "defs.h"
#include "param.h"
#include "buffer.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int  WriteReadIfValid      (int, int, uint64_t, uint64_t, uint64_t, int,
	              	   char **, char **, EBUF *, EBUF *, int64_t,
			   double, double);
void FilterCharacteristics (FC_PARAMETERS *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

