#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "defs.h"
#include "keys.h"
#include "msg.h"
#include "mem.h"
#include "dna.h"
#include "param.h"
#include "common.h"
#include "levels.h"
#include "strings.h"
#include "af.h"
#include "ma.h"
#include "fc.h"
#include "nc.h"
#include "ncd.h"
#include "rw.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NC_PARAMETERS  *NCP;
NCD_PARAMETERS *NCDP;
RW_PARAMETERS  *RWP;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_MapRaws(char **p, int c)
  {
  char    **xargv, *xpl = NULL;
  int32_t n, xargc = 0;

  RW_PARAMETERS *RWP = (RW_PARAMETERS *) Calloc(1, sizeof(RW_PARAMETERS));

  RWP->help     = ArgsState  (DEF_RW_HELP,    p, c, "-h", "--help");
  RWP->verbose  = ArgsState  (DEF_RW_VERBOSE, p, c, "-v", "--verbose");
  RWP->vv       = ArgsState  (DEF_RW_VV,      p, c, "-vv", "--very-verbose");
  RWP->force    = ArgsState  (DEF_RW_FORCE,   p, c, "-f", "--force");
  RWP->threads  = ArgsState  (DEF_RW_THREADS, p, c, "-t", "--threads");
  RWP->aa       = ArgsState  (DEF_RW_AA,      p, c, "-a", "--aminoacids");
  RWP->ir       = ArgsState  (DEF_RW_IR,      p, c, "-i", "--ignore-ir");
  RWP->out      = ArgsState  (DEF_RW_OUT,     p, c, "-o", "--stdout");
  RWP->plots    = ArgsState  (DEF_RW_PLOTS,   p, c, "-p", "--plots");
  RWP->min      = ArgsNum    (DEF_RW_MIN,     p, c, "-min", "--minimum",
                              DEF_RW_MIN_CTX, DEF_RW_MAX_CTX);
  RWP->max      = ArgsNum    (DEF_RW_MAX,     p, c, "-max", "--maximum",
                              DEF_RW_MIN_CTX, DEF_RW_MAX_CTX);

  if(c < MIN_NPARAM_FOR_PROGS+1 || RWP->help)
    {
    PrintMenuRW();
    return;
    }

  if(RWP->vv == 1) RWP->verbose = 1;

  CheckMinMax(RWP->min, RWP->max);

  RWP->nKmers = RWP->max - RWP->min + 1;

  if(RWP->threads != 0)
    RWP->threads = RWP->nKmers;

  RWP->parasite = CloneString(p[c-1]);
  CheckFileEmpty(RWP->parasite);

  RWP->host = CloneString(p[c-2]);
  CheckFileEmpty(RWP->host);

  RWP->prefix = CloneString(p[c-1]);

  if(RWP->verbose) PrintParametersRW(RWP);

  MapRAWs(RWP);

  if(RWP->verbose) fprintf(stderr, "[>] Done!\n");

  free(RWP);
  return;
  }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_NormalizedCompressionDistance(char **p, int c)
  {
  char    **xargv, *xpl = NULL;
  int32_t n, xargc = 0;

  NCDP = (NCD_PARAMETERS *) Calloc(1, sizeof(NCD_PARAMETERS));

  NCDP->help    = ArgsState  (DEF_NCD_HELP,     p, c, "-h", "--help");
  NCDP->verbose = ArgsState  (DEF_NCD_VERBOSE,  p, c, "-v", "--verbose");
  NCDP->dna     = ArgsState  (DEF_NCD_DNA,      p, c, "-d", "--dna");
  NCDP->threads = ArgsNum    (4,                p, c, "-t", "--threads",
                             DEF_NCD_MIN_THREADS, DEF_NCD_MAX_THREADS);
  NCDP->level   = ArgsNum    (0,                p, c, "-l", "--level",
                             DEF_NCD_MIN_LEVEL, DEF_NCD_MAX_LEVEL);

  if(c < MIN_NPARAM_FOR_PROGS + 1 || NCDP->help)
    {
    PrintMenuNCD();
    return;
    }

  if(ArgsState(0, p, c, "-s", "--show-levels"))
    {
    PrintLevelsNCD();
    exit(1);
    }

  if(ArgsState(0, p, c, "-p", "--show-parameters"))
    {
    PrintModels();
    exit(1);
    }

  NCDP->nModels = 0;
  for(n = 1 ; n < c ; ++n)
    if(strcmp(p[n], "-m") == 0)
      NCDP->nModels += 1;

  if(NCDP->nModels == 0 && NCDP->level == 0)
    NCDP->level = DEF_NCD_LEVEL;

  if(NCDP->level != 0)
    {
    xpl = GetLevelsLR(NCDP->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        NCDP->nModels += 1;
    }

  if(NCDP->nModels == 0)
    {
    PrintWarning("at least you need to use one model!");
    exit(1);
    }

  NCDP->model = (MODEL_PAR *) Calloc(NCDP->nModels, sizeof(MODEL_PAR));

  int k = 0;
  for(n = 1 ; n < c ; ++n)
    if(strcmp(p[n], "-m") == 0)
      NCDP->model[k++] = ArgsUniqModelNCD(p[n+1], 0);
  if(NCDP->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        NCDP->model[k++] = ArgsUniqModelNCD(xargv[n+1], 0);
    }

  k = 0;
  int ref_exists = 0;
  for(n = 1 ; n < c ; ++n)
    if(strcmp(p[n], "-r") == 0)
      {
      NCDP->reference = CloneString(p[n+1]);
      CheckFileIsFASTA(NCDP->reference);
      ref_exists = 1;
      break;
      }

  if(ref_exists == 0)
    {
    PrintWarning("reference does not exists!");
    exit(1);
    }

  NCDP->filename = p[c-1];
  CheckFileEmpty (NCDP->filename);

  if(NCDP->verbose) PrintParametersNCD(NCDP);

  NormalizedCompressionDistance(NCDP);

  if(NCDP->verbose) fprintf(stderr, "[>] Done!                  \n");

  free(NCDP);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_NormalizedCompression(char **p, int c)
  {
  char    **xargv, *xpl = NULL;
  int32_t n, xargc = 0;

  NCP = (NC_PARAMETERS *) Calloc(1, sizeof(NC_PARAMETERS));

  NCP->help    = ArgsState  (DEF_NC_HELP,     p, c, "-h", "--help");
  NCP->verbose = ArgsState  (DEF_NC_VERBOSE,  p, c, "-v", "--verbose");
  NCP->dna     = ArgsState  (DEF_NC_DNA,      p, c, "-d", "--dna");
  NCP->threads = ArgsNum    (4,               p, c, "-t", "--threads",
                            DEF_NC_MIN_THREADS, DEF_NC_MAX_THREADS);
  NCP->level   = ArgsNum    (0,               p, c, "-l", "--level",
                            DEF_NC_MIN_LEVEL, DEF_NC_MAX_LEVEL);

  if(c < MIN_NPARAM_FOR_PROGS + 1 || NCP->help)
    {
    PrintMenuNC();
    return;
    }

  if(ArgsState(0, p, c, "-s", "--show-levels"))
    {
    PrintLevelsNC();
    exit(1);
    }

  if(ArgsState(0, p, c, "-p", "--show-parameters"))
    {
    PrintModels();
    exit(1);
    }

  NCP->nModels = 0;
  for(n = 1 ; n < c ; ++n)
    if(strcmp(p[n], "-m") == 0)
      NCP->nModels += 1;

  if(NCP->nModels == 0 && NCP->level == 0)
    NCP->level = DEF_NC_LEVEL;

  if(NCP->level != 0)
    {
    xpl = GetLevelsLR(NCP->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        NCP->nModels += 1;
    }

  if(NCP->nModels == 0)
    {
    PrintWarning("at least you need to use a context model!");
    exit(1);
    }

  NCP->model = (MODEL_PAR *) Calloc(NCP->nModels, sizeof(MODEL_PAR));

  int k = 0;
  for(n = 1 ; n < c ; ++n)
    if(strcmp(p[n], "-m") == 0)
      NCP->model[k++] = ArgsUniqModelNC(p[n+1], 0);
  if(NCP->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        NCP->model[k++] = ArgsUniqModelNC(xargv[n+1], 0);
    }

  NCP->filename = p[c-1];
  CheckFileEmpty (NCP->filename);

  if(NCP->verbose) PrintParametersNC(NCP);

  NormalizedCompression(NCP);

  if(NCP->verbose) fprintf(stderr, "[>] Done!\n");

  free(NCP);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_FilterCharacteristics(char **p, int c)
  {
  FC_PARAMETERS *MAP = (FC_PARAMETERS *) Calloc(1, sizeof(FC_PARAMETERS));

  MAP->help     = ArgsState  (DEF_FC_HELP,     p, c, "-h",   "--help");
  MAP->verbose  = ArgsState  (DEF_FC_VERBOSE,  p, c, "-v",   "--verbose");
  MAP->complete = ArgsState  (DEF_FC_COMPLETE, p, c, "-c",   "--complete");
  MAP->minimum  = ArgsNum    (DEF_FC_MINIMUM,  p, c, "-min", "--minimum",
                             1, UINT_MAX);
  MAP->maximum  = ArgsNum    (UINT_MAX,        p, c, "-max", "--maximum",
                             1, UINT_MAX);
  MAP->alphabet = ArgsString (DEF_FC_ALPHABET, p, c, "-a",   "--alphabet");
  MAP->cg_min   = ArgsDouble (0.0,             p, c, "-ncg", "--cg-minimum");
  MAP->cg_max   = ArgsDouble (1.0,             p, c, "-mcg", "--cg-maximum");

  uint32_t x, y;

  MAP->nPatterns = 0;
  for(x = 2 ; x < c ; ++x)
    if(!strcmp(p[x], "-p") || !strcmp(p[x], "--pattern"))
      MAP->nPatterns++;

  MAP->nIgnore = 0;
  for(x = 2 ; x < c ; ++x)
    if(!strcmp(p[x], "-i") || !strcmp(p[x], "--ignore"))
      MAP->nIgnore++;

  MAP->patterns = (char **) Calloc(MAP->nPatterns, sizeof(char *));
  MAP->ignore   = (char **) Calloc(MAP->nIgnore,   sizeof(char *));

  y = 0;
  for(x = 2 ; x < c ; ++x)
    if(!strcmp(p[x], "-p") || !strcmp(p[x], "--pattern"))
      {
      int pLen = strlen(p[x+1]);
      MAP->patterns[y] = (char *) malloc(pLen + 1 * sizeof(char));
      strcpy(MAP->patterns[y], p[x+1]);
      MAP->patterns[y++][pLen + 1] = '\0';
      }

  y = 0;
  for(x = 2 ; x < c ; ++x)
    if(!strcmp(p[x], "-i") || !strcmp(p[x], "--ignore"))
      {
      int pLen = strlen(p[x+1]);
      MAP->ignore[y] = (char *) malloc(pLen + 1 * sizeof(char));
      strcpy(MAP->ignore[y], p[x+1]);
      MAP->ignore[y++][pLen + 1] = '\0';
      }

  if(c < MIN_NPARAM_FOR_PROGS || MAP->help)
    {
    PrintMenuFC();
    return;
    }

  CheckStdinEmpty();

  if(MAP->verbose) PrintParametersFC(MAP);

  FilterCharacteristics(MAP);

  if(MAP->verbose) fprintf(stderr, "[>] Done!\n");

  for(x = 0 ; x < MAP->nPatterns ; ++x)
    free(MAP->patterns[x]);
  free(MAP->patterns);

  for(x = 0 ; x < MAP->nIgnore ; ++x)
    free(MAP->ignore[x]);
  free(MAP->ignore);

  free(MAP);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_AlphabetFreqs(char **p, int c)
  {
  AF_PARAMETERS *MAP = (AF_PARAMETERS *) Calloc(1, sizeof(AF_PARAMETERS));

  MAP->help     = ArgsState  (DEF_AF_HELP,     p, c, "-h", "--help");
  MAP->verbose  = ArgsState  (DEF_AF_VERBOSE,  p, c, "-v", "--verbose");
  MAP->f_line   = ArgsState  (DEF_AF_FL,       p, c, "-p", "--first-line");
  MAP->alphabet = ArgsString (DEF_AF_ALPHABET, p, c, "-a", "--alphabet");

  if(c < MIN_NPARAM_FOR_PROGS || MAP->help)
    {
    PrintMenuAF();
    return;
    }	

  CheckStdinEmpty();
  
  if(MAP->verbose) PrintParametersAF(MAP);

  AlphabetFreqs(MAP);

  if(MAP->verbose) fprintf(stderr, "[>] Done!\n");

  free(MAP);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void P_MovingAvg(char **p, int c)
  {
  MA_PARAMETERS *MAP = (MA_PARAMETERS *) Calloc(1, sizeof(MA_PARAMETERS));

  MAP->help        = ArgsState (DEF_MA_HELP,     p, c, "-h", "--help");
  MAP->verbose     = ArgsState (DEF_MA_VERBOSE,  p, c, "-v", "--verbose");
  MAP->ignore      = ArgsState (DEF_MA_IGNORE,   p, c, "-i", "--ignore");
  MAP->show_pos    = ArgsState (DEF_MA_POSITION, p, c, "-p", "--position");
  MAP->column      = ArgsNum   (DEF_MA_COLUMN,   p, c, "-c", "--column", 1, 
		               UINT_MAX);
  MAP->window_size = ArgsNum   (DEF_MA_WINDOW,   p, c, "-w", "--window", 1, 
		               UINT_MAX);

  if(c < MIN_NPARAM_FOR_PROGS || MAP->help)
    {
    PrintMenuMA();
    return;
    }

  CheckStdinEmpty();

  if(MAP->verbose) PrintParametersMA(MAP);

  MovingAverage(MAP);

  if(MAP->verbose) fprintf(stderr, "[>] Done!\n");

  free(MAP);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[])
  {
  char **p = *&argv;

  if(ArgsState(0, p, argc, "-V", "--version"))
    {
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_HELP, p, argc, "-h", "--help") && argc == 2 || argc < 2)
    {
    PrintMenu();
    return EXIT_SUCCESS;
    }

  switch(KeyString(argv[1]))
    {
    case K1: PrintMenu();                                   break;
    case K2: P_MovingAvg                     (argv, argc);  break;
    case K3: P_FilterCharacteristics         (argv, argc);  break;
    case K4: P_AlphabetFreqs                 (argv, argc);  break;
    case K5: P_NormalizedCompression         (argv, argc);  break;
    case K6: P_NormalizedCompressionDistance (argv, argc);  break;
    case K7: P_MapRaws                       (argv, argc);  break;

    default:
    PrintWarning("unknown menu option!");
    PrintMenu();
    }

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
