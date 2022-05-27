#ifndef RW_H_INCLUDED
#define RW_H_INCLUDED

#include <stdio.h>
#include <stdint.h>

#include "defs.h"
#include "param.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void     CheckMinMax   (uint32_t, uint32_t);
void     MapRAWs       (RW_PARAMETERS *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
