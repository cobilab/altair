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
// - - - - - - A L P H A B E T   F R E Q U E N C I E S   M A I N - - - - - - -
//
void AlphabetFreqs(AF_PARAMETERS *MAP)
  {
  uint64_t x = 0;

  CheckInputIsFASTA();

  if(MAP->verbose) fprintf(stderr, "[>] Loading alphabet ...\n");
   
  char *symbols = (char *) Calloc(256, sizeof(char));  
  char *rev     = (char *) Calloc(256, sizeof(char));  

  int str_len = strlen(MAP->alphabet);
  if(str_len < 2)
    {
    PrintWarning("alphabet must contaim more than 1 symbol");
    exit(1);
    }

  StringHasUniqueCharacters(MAP->alphabet);

  for(x = 0 ; x < str_len ; ++x)
    {
    symbols[MAP->alphabet[x]] = 1;
    rev[MAP->alphabet[x]] = x;
    }

  if(MAP->verbose) 
    fprintf(stderr, "[>] Alphabet: ");

  unsigned nSymbols = 0;
  for(x = 0 ; x < 255 ; ++x)
    if(symbols[x] == 1)
      {
      ++nSymbols;
      if(MAP->verbose) fprintf(stderr, "%c", (int) x);
      }
  
  uint64_t *counts = (uint64_t *) Calloc(nSymbols + 1, sizeof(uint64_t));  
  
  if(MAP->verbose) fprintf(stderr, "\n");
  if(MAP->verbose) fprintf(stderr, "[>] Alphabet size: %u\n", nSymbols);
  if(MAP->verbose) fprintf(stderr, "[>] Running frequency ...\n");

  int64_t k, idx_in, seq_len = 0, nReads = 0;
  uint8_t *in_Buf, sym, header = 0;

  in_Buf = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  
  EBUF *header_Buf   = CreateEBuffer(5000,  200000);
  EBUF *sequence_Buf = CreateEBuffer(30000, 200000);

  nReads = 0;
  header = 1;
  seq_len = 0;

  if(MAP->f_line)
    {
    for(x = 0 ; x < nSymbols ; ++x)
      fprintf(stdout, "%c\t", MAP->alphabet[x]);
    fprintf(stdout, "\n");
    }

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
	  if(nReads > 0)
	    {
            for(x = 0 ; x < nSymbols ; ++x)
	      {
              fprintf(stdout, "%lf\t", (double) counts[x] / seq_len);
              counts[x] = 0;
	      }
            fprintf(stdout, "\n");
	    }

	  header = 0;
	  seq_len = 0;
	  }

	continue;
        }      
      else if(header == 0)
	{
	if(sym == '>')
          {
	  ++nReads;

	  ResetEBuffer(header_Buf);
	  ResetEBuffer(sequence_Buf);

	  header = 1;
	  goto anotherheader;
	  }
		
	sequence_Buf->buf[sequence_Buf->idx] = sym;
        UpdateEBuffer(sequence_Buf);

        if(symbols[sym] == 1)
	  {
	  ++seq_len;
          counts[rev[sym]]++;
	  }
	}
      else 
	{
        PrintWarning("exception found!");
	exit(1);
	}
      }

  ++nReads;
  for(x = 0 ; x < nSymbols ; ++x)
    fprintf(stdout, "%lf\t", (double) counts[x] / seq_len);
  fprintf(stdout, "\n");
  
  if(MAP->verbose) 
    fprintf(stderr, "[>] Number of processed FASTA reads: %"PRIi64"\n", nReads);
 
  Free(in_Buf);
  RemoveEBuffer(header_Buf); 
  RemoveEBuffer(sequence_Buf); 
  
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

