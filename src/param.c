#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>

#include "mem.h"
#include "common.h"
#include "param.h"
#include "msg.h"
#include "dna.h"
#include "strings.h"
#include "keys.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t ArgsNum(uint32_t d, char *a[], uint32_t n, char *s, char *s2,
uint32_t l, uint32_t u)
  {
  uint32_t x;
  for( ; --n ; ) if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
    {
    if((x = atol(a[n+1])) < l || x > u)
      {
      fprintf(stderr, "[x] Invalid number! Interval: [%u;%u].\n", l, u);
      exit(EXIT_FAILURE);
      }
    return x;
    }
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double ArgsDouble(double d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      return atof(a[n+1]);
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t ArgsState(uint8_t d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      return d == 0 ? 1 : 0;
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *ArgsString(char *d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      {
      if(a[n+1] == NULL || strlen(a[n+1]) < 1)
        {
        PrintWarning("string is empty!");
	exit(1);
        }
      return a[n+1];
      }
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MODEL_PAR ArgsUniqModelNC(char *str, uint8_t type)
  {
  uint32_t  ctx, den, ir, hash, edits, eDen, memory, models_exit = 0;
  double    gamma, eGamma;
  MODEL_PAR Mp;


  if(sscanf(str, "%u:%u:%u:%u:%u:%lf/%u:%u:%lf", &ctx, &den, &memory, &ir,
    &hash, &gamma, &edits, &eDen, &eGamma) == 9)
    {

    if(ctx > NC_MAX_CTX || ctx < NC_MIN_CTX)
      {
      fprintf(stderr, "ERROR: Invalid Context!\n");
      models_exit = 1;
      }

    if(den > NC_MAX_DEN || den < NC_MIN_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator!\n");
      models_exit = 1;
      }

    if(memory > NC_MAX_MEM)
      {
      fprintf(stderr, "ERROR: Invalid cache memory model!\n");
      models_exit = 1;
      }

    if(ir > 2)
      {
      fprintf(stderr, "ERROR: Invalid IR!\n");
      models_exit = 1;
      }

    if(hash > 255)
      {
      fprintf(stderr, "ERROR: Invalid cache-hash size!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma!\n");
      models_exit = 1;
      }

    if(edits > 20)
      {
      fprintf(stderr, "ERROR: Invalid number of editions (substitutions)!\n");
      models_exit = 1;
      }

    if(eDen > NC_MAX_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator (substitutions)!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma (substitutions)!\n");
      models_exit = 1;
      }

    if(models_exit == 1)
      {
      PrintModels();
      fprintf(stderr, "\nReset models according to the above description.\n");
      exit(1);
      }

    Mp.ctx      = ctx;
    Mp.den      = den;
    Mp.memory   = memory;
    Mp.ir       = ir;
    Mp.hashSize = hash;
    Mp.gamma    = ((int)( gamma * 65534)) / 65534.0;
    Mp.eGamma   = ((int)(eGamma * 65534)) / 65534.0;
    Mp.edits    = edits;
    Mp.eDen     = eDen;
    Mp.type     = type;
    return Mp;
    }
  else{
    fprintf(stderr, "Error: unknown scheme for model arguments!\n");
    PrintModels();
    fprintf(stderr, "\nReset models according to the above description.\n");
    exit(1);
    }

  return Mp;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MODEL_PAR ArgsUniqModelNCD(char *str, uint8_t type)
  {
  uint32_t  ctx, den, ir, hash, edits, eDen, memory, models_exit = 0;
  double    gamma, eGamma;
  MODEL_PAR Mp;


  if(sscanf(str, "%u:%u:%u:%u:%u:%lf/%u:%u:%lf", &ctx, &den, &memory, &ir,
    &hash, &gamma, &edits, &eDen, &eGamma) == 9)
    {

    if(ctx > NCD_MAX_CTX || ctx < NCD_MIN_CTX)
      {
      fprintf(stderr, "ERROR: Invalid Context!\n");
      models_exit = 1;
      }

    if(den > NCD_MAX_DEN || den < NCD_MIN_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator!\n");
      models_exit = 1;
      }

    if(memory > NCD_MAX_MEM)
      {
      fprintf(stderr, "ERROR: Invalid cache memory model!\n");
      models_exit = 1;
      }

    if(ir > 2)
      {
      fprintf(stderr, "ERROR: Invalid IR!\n");
      models_exit = 1;
      }

    if(hash > 255)
      {
      fprintf(stderr, "ERROR: Invalid cache-hash size!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma!\n");
      models_exit = 1;
      }

    if(edits > 20)
      {
      fprintf(stderr, "ERROR: Invalid number of editions (substitutions)!\n");
      models_exit = 1;
      }

    if(eDen > NCD_MAX_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator (substitutions)!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma (substitutions)!\n");
      models_exit = 1;
      }

    if(models_exit == 1)
      {
      PrintModels();
      fprintf(stderr, "\nReset models according to the above description.\n");
      exit(1);
      }

    Mp.ctx      = ctx;
    Mp.den      = den;
    Mp.memory   = memory;
    Mp.ir       = ir;
    Mp.hashSize = hash;
    Mp.gamma    = ((int)( gamma * 65534)) / 65534.0;
    Mp.eGamma   = ((int)(eGamma * 65534)) / 65534.0;
    Mp.edits    = edits;
    Mp.eDen     = eDen;
    Mp.type     = type;
    return Mp;
    }
  else{
    fprintf(stderr, "Error: unknown scheme for model arguments!\n");
    PrintModels();
    fprintf(stderr, "\nReset models according to the above description.\n");
    exit(1);
    }

  return Mp;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *ArgsFiles(char *arg[], uint32_t argc, char *str)
  {
  int32_t n = argc;

  for( ; --n ; )
    if(!strcmp(str, arg[n]))
      return CloneString(arg[n+1]);

  return Cat(Cat(arg[argc-2], arg[argc-1]), ".svg");
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersMA(MA_PARAMETERS *M)
  {
  fprintf(stderr, "[>] Running %s %s ...\n", PNAME, LT_KEYS[1].key);
  fprintf(stderr, "[>] Verbose mode: %s\n", M->verbose ? "yes": "no");
  fprintf(stderr, "[>] Ignore 1st line: %s\n", M->ignore ? "yes": "no");
  fprintf(stderr, "[>] Show positions: %s\n", M->show_pos ? "yes": "no");
  fprintf(stderr, "[>] Selected column: %u\n", M->column);
  fprintf(stderr, "[>] Window size: %u\n", M->window_size);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersAF(AF_PARAMETERS *M)
  {
  fprintf(stderr, "[>] Running %s %s ...\n", PNAME, LT_KEYS[3].key);
  fprintf(stderr, "[>] Verbose mode: %s\n", M->verbose ? "yes": "no");
  fprintf(stderr, "[>] Print 1st line: %s\n", M->f_line ? "yes": "no");
  fprintf(stderr, "[>] Alphabet: %s\n", M->alphabet);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersFC(FC_PARAMETERS *M)
  {
  fprintf(stderr, "[>] Running %s %s ...\n", PNAME, LT_KEYS[2].key);
  fprintf(stderr, "[>] Verbose mode: %s\n", M->verbose ? "yes": "no");
  fprintf(stderr, "[>] Complete alphabet: %s\n", M->complete ? "yes": "no");
  fprintf(stderr, "[>] Alphabet: %s\n", M->alphabet);
  if(M->nPatterns > 0) 
  fprintf(stderr, "[>] Considering patterns: \n");
  for(int x = 0 ; x < M->nPatterns ; ++x)
    fprintf(stderr, "[>] - %s\n", M->patterns[x]);
  if(M->nIgnore > 0) 
    fprintf(stderr, "[>] Considering ignore patterns: \n");
  for(int x = 0 ; x < M->nIgnore ; ++x)
    fprintf(stderr, "[>] - %s\n", M->ignore[x]);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersNC(NC_PARAMETERS *M)
  {
  uint32_t n;

  fprintf(stderr, "[>] Number of threads .............. %u\n", M->threads);
  fprintf(stderr, "[>] DNA mode ....................... %s\n", M->dna == 0 ?
  "no" : "yes");
  for(n = 0 ; n < M->nModels ; ++n){
    fprintf(stderr, "[>] Target model %d:\n", n+1);
    fprintf(stderr, "  [+] Context order ................ %u\n",
    M->model[n].ctx);
    fprintf(stderr, "  [+] Alpha denominator ............ %u\n",
    M->model[n].den);
    fprintf(stderr, "  [+] Cache memory model ........... %u\n",
    M->model[n].memory);
    switch(M->model[n].ir){
      case 0:
      fprintf(stderr, "  [+] Inverted repeats ............. no (regular)\n");
      break;
      case 1:
      fprintf(stderr, "  [+] Inverted repeats ............. mix (reg & inv)\n");
      break;
      case 2:
      fprintf(stderr, "  [+] Inverted repeats ............. inverted only\n");
      break;
      }
    fprintf(stderr, "  [+] Cache-hash size .............. %u\n",
    M->model[n].hashSize);
    fprintf(stderr, "  [+] Gamma ........................ %.3lf\n",
    M->model[n].gamma);
    fprintf(stderr, "  [+] Allowable substitutions ...... %u\n",
    M->model[n].edits);
    if(M->model[n].edits != 0){
      fprintf(stderr, "  [+] Substitutions alpha den ...... %u\n",
      M->model[n].eDen);
      fprintf(stderr, "  [+] Substitutions gamma .......... %.3lf\n",
      M->model[n].eGamma);
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersNCD(NCD_PARAMETERS *M)
  {
  uint32_t n;

  fprintf(stderr, "[>] Number of threads .............. %u\n", M->threads);
  fprintf(stderr, "[>] DNA mode ....................... %s\n", M->dna == 0 ?
  "no" : "yes");
  for(n = 0 ; n < M->nModels ; ++n){
    fprintf(stderr, "[>] Target model %d:\n", n+1);
    fprintf(stderr, "  [+] Context order ................ %u\n",
    M->model[n].ctx);
    fprintf(stderr, "  [+] Alpha denominator ............ %u\n",
    M->model[n].den);
    fprintf(stderr, "  [+] Cache memory model ........... %u\n",
    M->model[n].memory);
    switch(M->model[n].ir){
      case 0:
      fprintf(stderr, "  [+] Inverted repeats ............. no (regular)\n");
      break;
      case 1:
      fprintf(stderr, "  [+] Inverted repeats ............. mix (reg & inv)\n");
      break;
      case 2:
      fprintf(stderr, "  [+] Inverted repeats ............. inverted only\n");
      break;
      }
    fprintf(stderr, "  [+] Cache-hash size .............. %u\n",
    M->model[n].hashSize);
    fprintf(stderr, "  [+] Gamma ........................ %.3lf\n",
    M->model[n].gamma);
    fprintf(stderr, "  [+] Allowable substitutions ...... %u\n",
    M->model[n].edits);
    if(M->model[n].edits != 0){
      fprintf(stderr, "  [+] Substitutions alpha den ...... %u\n",
      M->model[n].eDen);
      fprintf(stderr, "  [+] Substitutions gamma .......... %.3lf\n",
      M->model[n].eGamma);
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintParametersRW(RW_PARAMETERS *M)
  {
  uint32_t n;
  fprintf(stderr, "[>] Running %s %s ...\n", PNAME, LT_KEYS[6].key);
  fprintf(stderr, "[>] Verbose mode: %s\n", M->verbose ? "yes" : "no");
  fprintf(stderr, "[>] Very verbose mode: %s\n", M->vv ? "yes" : "no");
  fprintf(stderr, "[>] Force mode: %s\n", M->force ? "yes" : "no");
  if(M->threads != 0)
    fprintf(stderr, "[>] Number of threads: %u\n", M->threads);
  fprintf(stderr, "[>] Aminoacids: %s\n", !M->aa ? "yes" : "no");
  fprintf(stderr, "[>] K-mer models:\n");
  for(n = 0 ; n < M->nKmers ; ++n)
    {
    fprintf(stderr, "    [>] K-mer model %u:\n", n+1);
    fprintf(stderr, "        [>] K-mer .................. %u\n", M->min+n);
    fprintf(stderr, "        [>] Use inversions ......... %s\n", !M->ir ?
    "no" : "yes");
    }
  fprintf(stderr, "[>] Plots: %s\n", M->plots ? "yes" : "no");
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

