#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <limits.h>
#include <ctype.h>

#include "defs.h"
#include "af.h"
#include "mem.h"
#include "msg.h"
#include "dna.h"
#include "buffer.h"
#include "common.h"
#include "strings.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - W R I T E   R E A D   I F   V A L I D - - - - - - - - - -
//
int WriteReadIfValid(int nPatterns, int nIgnore, uint64_t min, uint64_t max, 
  uint64_t sl, int ri, char **patterns, char **ignore, EBUF *HB, EBUF *SB, 
  int64_t cg_count, double cg_min, double cg_max)
  {
  unsigned rwrite = 1, x;

  if(nPatterns == 0) rwrite = 0;

  for(x = 0 ; x < nPatterns ; ++x)
    if(strcasestr((char *) HB->buf, (char *) patterns[x]))
      rwrite = 0;

  for(x = 0 ; x < nIgnore ; ++x)
    if(strcasestr((char *) HB->buf, (char *) ignore[x]))
      rwrite = 1;
 
  double cgv = (double) cg_count / sl; 
  if(sl > min && sl < max && cgv >= cg_min && cgv <= cg_max && ri == 0 
  && rwrite == 0) 
    {
    for(x = 0 ; x < HB->idx ; ++x) fprintf(stdout, "%c", HB->buf[x]);
    for(x = 0 ; x < SB->idx ; ++x) fprintf(stdout, "%c", SB->buf[x]);
   
    return 1;
    }
  return 0;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - F I L T E R   C H A R A C T E R I S T I C S - - - - - - - -
//
void FilterCharacteristics(FC_PARAMETERS *MAP)
  {
  uint64_t x = 0, y = 0, valid = 0;

  CheckInputIsFASTA();

  if(MAP->verbose) fprintf(stderr, "[>] Loading alphabet ...\n");
   
  char *symbols = (char *) Calloc(256, sizeof(char));  

  int str_len = strlen(MAP->alphabet);
  if(str_len < 2)
    {
    PrintWarning("alphabet must contain more than 1 symbol");
    exit(1);
    }
  
  StringHasUniqueCharacters(MAP->alphabet);

  for(x = 0 ; x < str_len ; ++x)
    symbols[MAP->alphabet[x]] = 1;

  if(MAP->verbose) fprintf(stderr, "[>] Alphabet: ");
  unsigned nSymbols = 0;
  for(x = 0 ; x < 255 ; ++x)
    if(symbols[x] == 1)
      {
      ++nSymbols;
      if(MAP->verbose) 
        fprintf(stderr, "%c", (int) x);
      }
  if(MAP->verbose) fprintf(stderr, "\n");

  if(MAP->verbose) 
    {	  
    fprintf(stderr, "[>] Alphabet size: %u\n", nSymbols);
    fprintf(stderr, "[>] Minimum sequence size: %u\n", MAP->minimum);
    fprintf(stderr, "[>] Maximum sequence size: %u\n", MAP->maximum);
    if(MAP->cg_min != 0.0)
      fprintf(stderr, "[>] Minimum CG: %lf\n", MAP->cg_min);
    if(MAP->cg_max != 1.0)
      fprintf(stderr, "[>] Maximum CG: %lf\n", MAP->cg_max);
    }

  if(MAP->verbose) fprintf(stderr, "[>] Running filter ...\n");

  int64_t k, idx_in, seq_len = 0, nReads = 0, cg_count = 0;
  uint8_t *in_Buf, sym, header = 0, skip = 0, ignore_read = 0, rwrite;

  in_Buf = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  
  EBUF *header_Buf = CreateEBuffer(5000, 200000);
  EBUF *sequence_Buf = CreateEBuffer(30000, 200000);

  cg_count = 0;
  nReads = 0;
  header = 1;
  seq_len = 0;
  while((k = fread(in_Buf, 1, BUFFER_SIZE, stdin)))
    for(idx_in = 0 ; idx_in < k ; ++idx_in)
      {
      sym = in_Buf[idx_in];
    
      if(header == 1)
	{
	anotherheader:
        header_Buf->buf[header_Buf->idx] = sym;
        UpdateEBuffer(header_Buf);

	if(sym == '\n')
	  {
	  header = 0;
	  seq_len = 0;
	  ignore_read = 0;
	  cg_count = 0;
	  }

	continue;
        }      
      else if(header == 0)
	{
	if(sym == '>')
          {
	  ++nReads;

          if(WriteReadIfValid(MAP->nPatterns, MAP->nIgnore, MAP->minimum, 
	  MAP->maximum, seq_len, ignore_read, MAP->patterns, MAP->ignore, 
	  header_Buf, sequence_Buf, cg_count, MAP->cg_min, MAP->cg_max) == 1)
            ++valid; 

	  ResetEBuffer(header_Buf);
	  ResetEBuffer(sequence_Buf);

	  header = 1;
	  goto anotherheader;
	  }
		
	sequence_Buf->buf[sequence_Buf->idx] = sym;
        UpdateEBuffer(sequence_Buf);

        if(symbols[sym] == 1) ++seq_len;
        if(symbols[sym] == 0 && MAP->complete == 1) ignore_read = 1;
	if(sym == 'C' || sym == 'G')
	  ++cg_count;
	}
      else 
	{
        PrintWarning("exception found!");
	exit(1);	
	}
      }

  ++nReads;
  if(WriteReadIfValid(MAP->nPatterns, MAP->nIgnore, MAP->minimum,
  MAP->maximum, seq_len, ignore_read, MAP->patterns, MAP->ignore,
  header_Buf, sequence_Buf, cg_count, MAP->cg_min, MAP->cg_max) == 1)
    ++valid;
  
  if(MAP->verbose)
    {
    fprintf(stderr, "[>] Number of FASTA reads: %"PRIi64"\n", nReads);
    fprintf(stderr, "[>] Number of valid reads: %"PRIi64"\n", valid);
    }
 
  Free(in_Buf);
  RemoveEBuffer(header_Buf); 
  RemoveEBuffer(sequence_Buf); 
  
  if(MAP->verbose) fprintf(stderr, "[>] Done\n");

  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

