#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <pthread.h>

#include "defs.h"
#include "mem.h"
#include "msg.h"
#include "rw.h"

#include "param.h"
#include "common.h"
#include "buffer.h"
#include "dna.h"
#include "alphabet.h"
#include "kmer.h"
#include "dist.h"
#include "stats.h"
#include "strings.h"

RW_PARAMETERS *P;  // SHARE PARAMETERS FOR THREADS
ALPHABET      *AL; // SHARE PARAMETERS FOR THREADS
FASTA_READS   *FA; // SHARE PARAMETERS FOR THREADS

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - G E N E R A T E   G N U P L O T   S H - - - - - - - - -
//
/*
void PrintGnuplot(uint32_t min, uint32_t max, uint32_t max_nsym){

  uint32_t x;

  FILE *OUT = Fopen("kplot.sh", "w");
  char *colors[] = { "red", "royalblue", "grey", "green", "blue", 
	             "purple", "red", "yellow", "violet", "cyan", "azure", 
		     "chartreuse", "darkgreen", "firebrick", "goldenrod", 
		     "mediumaquamarine", "olive", "orangered", "navy", 
		     "orchid", "peru" };

  fprintf(OUT, "#!/bin/sh\n");
  fprintf(OUT, "echo 'set mapping cartesian\n");
  fprintf(OUT, "set output \"%s-k-map.pdf\"\n", P->prefix);
  fprintf(OUT, "set view 360,0,1,1 #0,0,1,1\n");
  fprintf(OUT, "set auto\n");
  fprintf(OUT, "set tics nomirror out scale 0.75\n");
  fprintf(OUT, "set zrange [%u:%u]\n", max, min);
  fprintf(OUT, "set xrange [0:%u]\n", max_nsym + 1);
  fprintf(OUT, "set yrange [0:%"PRIu64"]\n", FA->nReads+1);
  fprintf(OUT, "set ztics 1\n");
  fprintf(OUT, "set isosamples 60\n");
  fprintf(OUT, "set hidden3d\n");
  fprintf(OUT, "unset key\n");
  fprintf(OUT, "set palette defined (");
  for(x = min ; x <= max ; ++x){ 
    fprintf(OUT, "%u \"%s\"", x, colors[x-min]);
    if(x != max)
      fprintf(OUT, ", ");
    }
  fprintf(OUT, ")\n");
  fprintf(OUT, "set zlabel \"K-mer\"\n");
  fprintf(OUT, "set ylabel \"Strain\"\n");
  fprintf(OUT, "set xlabel \"Length\"\n");
  fprintf(OUT, "splot \"%s.mink\" u 2:1:3 with points pt 5 ps 0.2 palette'",
  P->prefix);
  fprintf(OUT, " | gnuplot -persist\n");

  return;
  }
*/
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - P R I N T   K   H E A D - - - - - - - - - - - - -
//

void PrintKHead(uint32_t k)
  {
  if(P->verbose)
    fprintf(stderr, "[>] K-mer %2u ...                       \n", k);
  return;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - A V E R A G E   C G   P R O F I L E - - - - - - - - -
//
/*
void AverageCGPlot(void)
  {
  uint32_t x, size = FA->max_nsym + 1;
  double v, *window = (double *) Calloc(size, sizeof(double));

  for(x = 0 ; x < P->nTar ; ++x)
     {
    char *name = (char *) Calloc(2048, sizeof(char));
    sprintf(name, "%s.cg.prof", P->tar[x]);
    FILE *IN = Fopen(name, "r");

    uint32_t idx = 0;
    while(fscanf(IN, "%lf\n", &v) != EOF){
      window[idx] += (double) v;
      if(++idx == size){
	fprintf(stderr, "[x] Error: out of bounds [AverageCGPlot]!\n");
	exit(1);
        }
      }
    Free(name);
    fclose(IN);
    }

  char *name = (char *) Calloc(2048, sizeof(char));
  sprintf(name, "%s.cg.avg", P->prefix);

  FILE *OUT = Fopen(name, "w");
  for(x = 0 ; x < size ; ++x)
    fprintf(OUT, "%.3lf\n", window[x] / FA->nReads);
  fclose(OUT);

  FILE *IN = Fopen(name, "r");
  char *name_out = (char *) Calloc(2048, sizeof(char));
  sprintf(name_out, "%s.cg.fil", P->prefix);
  FILE *OUT2 = Fopen(name_out, "w");
  Filter(IN, OUT2);
  fclose(OUT2);
  fclose(IN);

  if(P->plots)
    {
    FILE *OUT3 = Fopen("CGprofileplot.sh",  "w");

    fprintf(OUT3, "#!/bin/sh\n");
    fprintf(OUT3, "echo 'reset\n");
    fprintf(OUT3, "set terminal pdfcairo enhanced font \"Verdana,12\"\n");
    fprintf(OUT3, "set output \"CGprofile.pdf\"\n");
    fprintf(OUT3, "set style line 11 lc rgb \"#808080\" lt 1\n");
    fprintf(OUT3, "set border 3 back ls 11\n");
    fprintf(OUT3, "set tics nomirror\n");
    fprintf(OUT3, "set size ratio 0.1\n");
    fprintf(OUT3, "unset grid\n");
    fprintf(OUT3, "set style line 1 lc rgb \"#0060ad\" lt 1 lw 1\n");
    fprintf(OUT3, "unset key\n");
    fprintf(OUT3, "set xlabel \"Length\"\n");
    fprintf(OUT3, "set ylabel \"CG %%\"\n");
    fprintf(OUT3, "set xrange [%u:%u]\n", 0, size);
    fprintf(OUT3, "set yrange [%u:%u]\n", 0, 1);
    fprintf(OUT3, "plot \"%s\" u 1 w lines ls 1' | gnuplot -persist\n", 
		  name_out);

    fclose(OUT3);
    if(P->verbose)
      fprintf(stderr, "[>] Shell plot with CG profile %% in CGprofileplot.sh\n");
    }

  Free(name_out);
  Free(name);

  return;
  }
*/
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - P R I N T   C G   P R O F I L E - - - - - - - - - -
//

void PrintCGProfile(uint32_t id){

  int v, size = 20;
  uint32_t x, idx = 0, bases = 0;
  uint32_t *window = (uint32_t *) Calloc(size +1, sizeof(uint32_t));

  char *name = (char *) Calloc(2048, sizeof(char));
  sprintf(name, "%s.cg.prof", P->prefix);
  FILE *OUT = Fopen(name, "w"); 
  FILE *IN  = Fopen(P->prefix, "r"); 

  while((v = fgetc(IN)) != EOF){
    window[idx] = v;
    uint32_t ALL = 0, AT = 0, CG = 0;
    for(x = 0 ; x < size ; ++x){
      switch(window[x]){
	case 'A': case 'a': case 'T': case 't': AT++; ALL++; break;
	case 'C': case 'c': case 'G': case 'g': CG++; ALL++; break;
	default: break;
        }
      }
    if(++bases > size){
      if(AT == 0 || CG == 0 || ALL == 0)
        fprintf(OUT, "%.4lf\n", (double) 0.0);  
      else	      
        fprintf(OUT, "%.4lf\n", (double) CG / ALL);  
      } 
    if(++idx == size)
      idx = 0;    
    }

  fclose(IN);
  fclose(OUT);
  Free(window);
  Free(name);

  return;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - C R E A T E   M O D E L S - - - - - - - - - - - -
//

KModels *CreateKModels(uint32_t nModels)
  {
  uint32_t x;
  KModels *KM = (KModels *) Calloc(1, sizeof(KModels));
  KM->nModels = nModels;
  KM->M       = (KModel **) Calloc(KM->nModels, sizeof(KModel *));
  for(x = 0 ; x < KM->nModels ; ++x)
    {
    KM->M[x]  = CreateKModel(P->min + x, P->ir, P->aa, P->nSym);
    KM->M[x]->id = x;
    }
  return KM;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - M I N   M A X   C H E C K - - - - - - - - - - - -
//

void CheckMinMax(uint32_t min, uint32_t max)
  {
  if(min > max)
    {
    fprintf(stderr, "[x] Error: minimum k-mer is higher than maximum!\n");
    exit(EXIT_FAILURE);
    }
  return;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - P R I N T   A V E R A G E   P L O T - - - - - - - - - - 
//

void PrintAveragePlot(Stats *S)
  {
  uint32_t x;

  FILE *OUT = Fopen("AVG-data.eg", "w");
  for(x = 0 ; x < P->nKmers ; ++x)
    fprintf(OUT, "%u\t%.3lf\n", P->min + x, S->average[x]);
  fclose(OUT);

  FILE *OUT2 = Fopen("AVGplot.sh", "w");
  fprintf(OUT2, "#!/bin/sh\n");
  fprintf(OUT2, "echo 'reset\n");
  fprintf(OUT2, "set terminal pdfcairo enhanced font \"Verdana,12\"\n");
  fprintf(OUT2, "set output \"AVGplot.pdf\"\n");
  //fprintf(OUT2, "set style line 11 lc rgb \"#808080\" lt 1\n");
  fprintf(OUT2, "set style line 11 lc rgb \"#000000\" lt 1\n");
  fprintf(OUT2, "set border 3 back ls 11\n");
  fprintf(OUT2, "set logscale y\n");
  fprintf(OUT2, "set tics nomirror\n");
  fprintf(OUT2, "set size ratio 1\n");
  //fprintf(OUT2, "set style line 12 lc rgb \"#808080\" lt 0 lw 2\n");
  fprintf(OUT2, "set style line 12 lc rgb \"#000000\" lt 0 lw 2\n");
  fprintf(OUT2, "set grid back ls 12\n");
  fprintf(OUT2, "set style line 1 lc rgb \"#8b1a0e\" pt 1 ps 1 lt 1 lw 3\n");
  fprintf(OUT2, "set style line 2 lc rgb \"#009900\" pt 7 ps 0.6 lt 1 lw 2\n");
  fprintf(OUT2, "unset key\n");
  fprintf(OUT2, "set xlabel \"k-mer\"\n");
  fprintf(OUT2, "set ylabel \"Number of mRAWs\"\n");
  fprintf(OUT2, "set xtics 1\n");
  fprintf(OUT2, "set xrange [%u:%u]\n", P->min, P->min + P->nKmers - 1);
  fprintf(OUT2, "set yrange [:]\n");
  fprintf(OUT2, "plot \"AVG-data.eg\" u 1:2 w lp ls 2' | gnuplot -persist\n");
  if(P->verbose) 
    fprintf(stderr, "[>] Shell plot with AVG %% in AVGplot.sh\n");
  fclose(OUT2);

  return;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - W R I T E   W O R D - - - - - - - - - - - - -
//

void RWord(FILE *F, uint8_t *b, int32_t i, uint32_t ctx)
  {
  uint8_t w[ctx+1], n;
  i -= ctx;
  for(n = 0 ; n < ctx ; ++n) w[n] = NumToDNASym(b[i+n]);
  w[ctx] = '\0';
  fprintf(F, "%s\n", w);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - W R I T E   W O R D   A A - - - - - - - - - - - -
//

void RWordAA(FILE *F, uint8_t *b, int32_t i, uint32_t ctx)
  {
  uint8_t w[ctx+1], n;
  i -= ctx;
  for(n = 0 ; n < ctx ; ++n)
    w[n] = AL->toChars[b[i+n]];
  w[ctx] = '\0';
  fprintf(F, "%s\n", w);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - T A R G E T   A A - - - - - - - - - - - - - -
//
// It finds and localizes the RAWS for a certain k and returns the number
// for stats print. Also, counts the number of symbols in each word. For AA.

uint32_t TargetAA(KModel *M, Dist *D, uint64_t id_tar, char *SEQ, uint32_t sz,
  FILE *Pos)
  {
  uint64_t i = 0, raw = 0, unknown = 0;
  uint64_t n, k, idxPos, hIndex;
  int64_t  idx = 0;
  uint8_t  *wBuf, sym, found = 0;
  CBUF     *sBuf;

  if(P->vv)
    fprintf(stderr, "[>] Searching target sequence %"PRIu64" ...\n", id_tar + 1);

  wBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  sBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);

  int raw_exists = 0;
  for(idxPos = 0 ; idxPos < sz ; ++idxPos)
    {
    ++i;
    sym = AL->revMap[SEQ[idxPos]];
    sBuf->buf[sBuf->idx] = sym;
    GetKIdxAA(sBuf->buf+sBuf->idx-1, M);

    if(i > M->ctx)
      {  // SKIP INITIAL CONTEXT, ALL "AAA..."
      if(M->mode == 0)
        { // TABLE MODE
        if(!M->array.counters[M->idx])
          { // THERE IS A RAW!
          fprintf(Pos, "%"PRIu64"\t%"PRIu64"\t", id_tar+1, i-M->ctx);
          RWordAA(Pos,sBuf->buf, sBuf->idx, M->ctx);
          raw_exists = 1;
          ++raw;
          }
        else
          { // NO RAW!
          raw_exists = 0;
          }
        }
      else
        { // HASH TABLE
        found = 0;
        hIndex = M->idx % K_HASH_SIZE;
        for(n = 0 ; n < M->hash.entrySize[hIndex] ; n++)
          if(((uint64_t) M->hash.keys[hIndex][n]*K_HASH_SIZE)+hIndex == M->idx)
            {
            found = 1;
            break;
            }
        if(found == 0)
          { // THERE IS A RAW
          fprintf(Pos, "%"PRIu64"\t%"PRIu64"\t", id_tar+1, i-M->ctx);
          RWordAA(Pos, sBuf->buf, sBuf->idx, M->ctx);
          raw_exists = 1;
          ++raw;
          }
        else
          { // NO RAW!
          raw_exists = 0;
          }
        }

      if(raw_exists == 1) // COUNTS THE BASE DISTRIBUTION OF THE RAWS
        for(n = 0 ; n < M->ctx ; n++)
          D->C[M->id].sym[sBuf->buf[sBuf->idx-n]] += 1;

      }

    UpdateCBuffer(sBuf);
    }

  ResetKIdx(M);
  Free(wBuf);
  RemoveCBuffer(sBuf);

  if(P->vv)
    fprintf(stdout, "[>] Kmer %u , Read %"PRIu64" , mRAWs: %.4lf %% ( %"PRIu64" in "
    "%"PRIu64" , unknown: %"PRIu64" , total: %u )\n", M->ctx, id_tar+1,
    (double) raw/(sz-unknown)*100.0, raw, sz-unknown, unknown, sz);

  if(P->vv)
    fprintf(stderr, "[>] Done!                           \n"); // VALID SPACES

  return raw;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - T A R G E T - - - - - - - - - - - - - - -
//
// It finds and localizes the RAWS for a certain k and returns the number
// for stats print. Also, counts the number of CG bases in each word.

uint32_t Target(KModel *M, Dist *D, uint64_t id_tar, char *SEQ, uint32_t sz,
  FILE *Pos)
  {
  uint64_t i = 0, raw = 0, unknown = 0;
  uint64_t n, k, idxPos, hIndex;
  int64_t  idx = 0;
  uint8_t  *wBuf, sym, found = 0;
  CBUF     *sBuf;

  if(P->vv)
    fprintf(stderr, "[>] Searching target sequence %"PRIu64" ...\n", id_tar+1);

  wBuf  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  sBuf  = CreateCBuffer(BUFFER_SIZE, BGUARD);

  int raw_exists = 0;
  for(idxPos = 0 ; idxPos < sz ; ++idxPos)
    {
    ++i;
    if((sym = DNASymToNum(SEQ[idxPos])) == 4)
      {
      ++unknown;
      continue;
      }
      
    sBuf->buf[sBuf->idx] = sym;
    GetKIdx(sBuf->buf+sBuf->idx-1, M); 
      
    if(i > M->ctx)
      {  // SKIP INITIAL CONTEXT, ALL "AAA..."
      if(M->mode == 0)
        { // TABLE MODE
        if(!M->array.counters[M->idx])
	  { // THERE IS A RAW!
          fprintf(Pos, "%"PRIu64"\t%"PRIu64"\t", id_tar+1,i-M->ctx);
          RWord(Pos,sBuf->buf, sBuf->idx, M->ctx);
	  raw_exists = 1;
          ++raw;
          }
        else
	  { // NO RAW!
	  raw_exists = 0;
          }
        }
      else
        { // HASH TABLE
        found = 0;
        hIndex = M->idx % K_HASH_SIZE;
        for(n = 0 ; n < M->hash.entrySize[hIndex] ; n++)
          if(((uint64_t) M->hash.keys[hIndex][n]*K_HASH_SIZE)+hIndex == M->idx)
	    {
            found = 1;
            break;
            }
        if(found == 0)
	  { // THERE IS A RAW
          fprintf(Pos, "%"PRIu64"\t%"PRIu64"\t", id_tar+1, i-M->ctx); 
          RWord(Pos, sBuf->buf, sBuf->idx, M->ctx);
	  raw_exists = 1;
          ++raw;
          }
        else
	  { // NO RAW!
	  raw_exists = 0;
          }
        }

      if(raw_exists == 1) // COUNTS THE BASE DISTRIBUTION OF THE RAWS
        for(n = 0 ; n < M->ctx ; n++)
	  D->C[M->id].sym[sBuf->buf[sBuf->idx-n]] += 1;

      }

    UpdateCBuffer(sBuf);
    }

  ResetKIdx(M);
  Free(wBuf);
  RemoveCBuffer(sBuf);

  if(P->vv)
    fprintf(stdout, "[>] Kmer %u , Read %"PRIu64" , mRAWs: %.4lf %% ( %"PRIu64""
    " in %"PRIu64" , unknown: %"PRIu64" , total: %u )\n", M->ctx, id_tar+1, 
    (double) raw/(sz-unknown)*100.0, raw, sz-unknown, unknown, sz);

  if(P->vv)
    fprintf(stderr, "[>] Done!                           \n"); // VALID SPACES

  return raw;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - R E F E R E N C E   A A - - - - - - - - - - - -
// 

void LoadReferenceAA(KModel *M){

  FILE     *Reader = Fopen(P->host, "r");
  int64_t  k, idxPos, header = 0;
  int64_t  idx = 0;
  uint8_t  *rBuf, *sBuf, sym;
  uint64_t i = 0;
  #ifdef PROGRESS
  uint64_t size = NBytesInFile(Reader);
  #endif

  if(P->verbose)
    {
    if(P->threads == 0)
      fprintf(stderr, "[>] Building host k-mer model (k=%u) ...\n", M->ctx);
    else if(P->threads == M->id+1)
      fprintf(stderr, "[>] Building host k-mer models ...\n");
    }

  rBuf  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  sBuf  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD+1, sizeof(uint8_t));
  sBuf += BGUARD;

  uint64_t skip = 0;
  while((k = fread(rBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      ++i;
      #ifdef PROGRESS
      if(P->threads == 0)              CalcProgress(size, i);
      else if(P->threads == M->id + 1) CalcProgress(size, i);
      #endif
      switch(rBuf[idxPos])
        {
        case '>':  header = 1; skip = 0; continue;
        case '\n': header = 0; continue;
        default: if(header==1) continue;
        }

      sym = AL->revMap[rBuf[idxPos]];
      sBuf[idx] = sym;
      GetKIdxAA(sBuf+idx-1, M);

      if(++skip > M->ctx)
        UpdateK(M); // SKIP INITIAL CONTEXT, ALL "AAA..."

      if(++idx == BUFFER_SIZE)
        {
        memcpy(sBuf-BGUARD, sBuf+idx-BGUARD, BGUARD);
        idx = 0;
        }
      }

  ResetKIdx(M);
  Free(rBuf);
  Free(sBuf-BGUARD);
  fclose(Reader);

  if(P->verbose)
    {
    if(P->threads == 0)
      fprintf(stderr, "[>] Done!                    \n");  // SPACES ARE VALID  
    else if(P->threads == M->id+1)
      fprintf(stderr, "[>] Done!                    \n");  // SPACES ARE VALID  
    }
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -
// 

void LoadReference(KModel *M)
  {
  FILE     *Reader = Fopen(P->host, "r");
  int64_t  k, idxPos, header = 0;
  int64_t  idx = 0;
  uint8_t  *rBuf, *sBuf, sym;
  uint64_t i = 0;
  #ifdef PROGRESS
  uint64_t size = NBytesInFile(Reader);
  #endif

  if(P->verbose)
    {
    if(P->threads == 0)
      fprintf(stderr, "[>] Building host k-mer model (k=%u) ...\n", M->ctx);
    else if(P->threads == M->id+1)
      fprintf(stderr, "[>] Building host k-mer models ...\n");
    }

  rBuf  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  sBuf  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD+1, sizeof(uint8_t));
  sBuf += BGUARD;

  uint64_t skip = 0;
  while((k = fread(rBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      ++i;
      #ifdef PROGRESS
      if(P->threads == 0)              CalcProgress(size, i);
      else if(P->threads == M->id + 1) CalcProgress(size, i);
      #endif
      switch(rBuf[idxPos])
        {
        case '>':  header = 1; skip = 0; continue;
        case '\n': header = 0; continue;  
        default: if(header==1) continue;
        }

      if((sym = DNASymToNum(rBuf[idxPos])) == 4) continue;    
      sBuf[idx] = sym;
      GetKIdx(sBuf+idx-1, M);

      if(++skip > M->ctx)
        { // SKIP INITIAL CONTEXT, ALL "AAA..."
        UpdateK(M);
        if(M->ir == 1)
	  {  // Inverted repeats
          GetKIdxIR(sBuf+idx, M);
          UpdateKIR(M);
          }
        }

      if(++idx == BUFFER_SIZE)
        {
        memcpy(sBuf-BGUARD, sBuf+idx-BGUARD, BGUARD);
        idx = 0;
        }
      }
 
  ResetKIdx(M);
  Free(rBuf);
  Free(sBuf-BGUARD);
  fclose(Reader);

  if(P->verbose)
    {
    if(P->threads == 0)
      fprintf(stderr, "[>] Done!                    \n"); // SPACES ARE VALID  
    else if(P->threads == M->id+1)
      fprintf(stderr, "[>] Done!                    \n"); // SPACES ARE VALID  
    }

  return;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - T H R E A D I N G - - - - - - - - - - - - - - -
//

void *LoadRefThread(void *Par)
  {
  KModel *M = (KModel *) Par;
  if(P->aa) LoadReference(M);
  else      LoadReferenceAA(M);
  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - E A G L E   M A I N - - - - - - - - - - - - - - -
//

void MapRAWs(RW_PARAMETERS *RWP){

  P = *&RWP; // THREADS SHARING
  uint32_t n, x, idx_model;
  uint64_t idx_reads;

  CheckFileIsFASTA(P->parasite);

  if(P->verbose) fprintf(stderr, "[>] Loading alphabet ...\n");
  AL = CreateAlphabet();
  LoadAlphabet2Files(AL, P->parasite, P->host);
  if(P->verbose)
    PrintAlphabet(AL);
  if(P->verbose) fprintf(stderr, "[>] Done!\n");

  if(P->aa){
    P->nSym = 4;
    if(AL->cardinality > 10){
      fprintf(stderr, "Warning: input files with large alphabet cardinality!\n");
      fprintf(stderr, "Tip: use -a flag to consider aminoacids.\n");
      fprintf(stderr, "Tip: use -f flag to force.\n");
      if(!P->force)
        exit(1);
      }
    }
  else{
    P->nSym = AL->cardinality;
    }

  if(P->verbose) fprintf(stderr, "[>] Analyzing parasite reads...\n");
  FA = CreateFastaReads();
  CountFastaReadsAndMax(FA, P->parasite);
  if(P->verbose)
    {
    fprintf(stderr, "[>] Number of parasite reads: %"PRIu64"\n", FA->nReads);
    fprintf(stderr, "[>] Maximum parasite length: %u\n", FA->max_nsym);
    }
  if(P->verbose) fprintf(stderr, "[>] Done!\n");
  
  if(P->verbose) fprintf(stderr, "[>] Initializing k-mer models ...\n");
  KModels *KM = CreateKModels(P->nKmers);
  if(P->verbose) fprintf(stderr, "[>] Done!\n");

  Stats    *S = CreateStats (P->nKmers);
  Dist     *D = CreateDist  (P->nKmers, P->nSym);

  char *SEQ = (char *) Calloc(FA->max_nsym + 1, sizeof(char));
  
  if(P->threads == 0){ // IF NO THREADS ==========================================

    if(P->verbose) fprintf(stderr, "[>] Using single-thread ...\n");

    for(idx_model = 0 ; idx_model < P->nKmers ; ++idx_model)
      {
      if(P->aa) LoadReference   (KM->M[idx_model]);
      else      LoadReferenceAA (KM->M[idx_model]);
     
      PrintKHead(idx_model + P->min);

      double sum = 0;
      uint64_t vr[FA->nReads+1];
      if(P->verbose) fprintf(stderr, "[>] Searching parasite sequences ...\n");

      char *name = (char *) Calloc(5048, sizeof(char));
      sprintf(name, "-k%u.eg", KM->M[idx_model]->ctx);
      char *name2 = Cat(P->prefix, name);
      FILE *Pos = Fopen(name2, "w");
      FILE *F = Fopen(P->parasite, "r");
      uint8_t buffer[BUFFER_SIZE], sym = 0, header = 1;
      uint32_t k, idx;
      uint64_t nSymbols = 0;

      idx_reads = 0;
      while((k = fread(buffer, 1, BUFFER_SIZE, F)))
        for(idx = 0 ; idx < k ; ++idx)
          {
          sym = buffer[idx];
          if(sym == '>')
	    { 
	    header = 1; 
	    if(idx_reads > 0)
	      {
              if(P->aa)
	        { 
	        vr[idx_reads-1] = Target(KM->M[idx_model], D, idx_reads-1, 
		                  SEQ, nSymbols, Pos);
	        }
              else
		{
                vr[idx_reads-1] = TargetAA(KM->M[idx_model], D, idx_reads-1, 
	                          SEQ, nSymbols, Pos);
		}
              sum += vr[idx_reads-1];
	      }
	    ++idx_reads; 
	    continue; 
	    }
          if(sym == '\n' && header == 1)
	    { header = 0; nSymbols = 0; continue; }
          if(sym == '\n') continue;
          if(header == 1) continue;
          SEQ[nSymbols++] = sym;
          }
      if(P->aa)
        vr[idx_reads-1] = Target(KM->M[idx_model], D, idx_reads-1,
                          SEQ, nSymbols, Pos);
      else
        vr[idx_reads-1] = TargetAA(KM->M[idx_model], D, idx_reads-1,
                          SEQ, nSymbols, Pos);

      fclose(F);
      fclose(Pos);
      Free(name);
      Free(name2);

      S->average[idx_model] = sum / (double) FA->nReads;

      if(P->verbose) fprintf(stderr, "[>] Done!                           "
      "                      \n"); // SPACES ARE VALID!

      double sum1 = 0;
      for(idx_reads = 0 ; idx_reads < FA->nReads ; ++idx_reads)
        sum1 = sum1 + pow(((double) vr[idx_reads] - S->average[idx_model]), 2);

      S->variance[idx_model] = sum1 / (double) FA->nReads;
      S->stddev  [idx_model] = sqrt(S->variance[idx_model]);

      DeleteKModel(KM->M[idx_model]);
      }

    if(P->plots) PrintAveragePlot(S);
    }
  else
    { // IF THREADS ===========================================================

    if(P->verbose) 
      fprintf(stderr, "[>] Using multi-threading to load host ...\n");

    pthread_t t[P->threads];
 
    for(n = 0 ; n < P->threads ; ++n)
      pthread_create(&(t[n+1]), NULL, LoadRefThread, (void *) KM->M[n]);
    for(n = 0 ; n < P->threads ; ++n) // DO NOT JOIN THIS FOR WITH PREVIOUS!
      pthread_join(t[n+1], NULL);
  
    for(idx_model = 0 ; idx_model < P->nKmers ; ++idx_model)
      {
      PrintKHead(idx_model + P->min);

      double sum = 0;
      uint64_t vr[FA->nReads+1];
      if(P->verbose) fprintf(stderr, "[>] Searching parasite sequences ...\n");

      char *name = (char *) Calloc(5048, sizeof(char));
      sprintf(name, "-k%u.eg", KM->M[idx_model]->ctx);
      char *name2 = Cat(P->prefix, name);
      FILE *Pos = Fopen(name2, "w");
      FILE *F = Fopen(P->parasite, "r");
      uint8_t buffer[BUFFER_SIZE], sym = 0, header = 1;
      uint32_t k, idx;
      uint64_t nSymbols = 0;

      idx_reads = 0;
      while((k = fread(buffer, 1, BUFFER_SIZE, F)))
        for(idx = 0 ; idx < k ; ++idx)
          {
          sym = buffer[idx];
          if(sym == '>')
            {
            header = 1;
            if(idx_reads > 0)
              {
              if(P->aa)
                {
                vr[idx_reads-1] = Target(KM->M[idx_model], D, idx_reads-1,
                                  SEQ, nSymbols, Pos);
                }
              else
                {
                vr[idx_reads-1] = TargetAA(KM->M[idx_model], D, idx_reads-1,
                                  SEQ, nSymbols, Pos);
                }
              sum += vr[idx_reads-1];
              }
            ++idx_reads;
            continue;
            }
          if(sym == '\n' && header == 1)
            { header = 0; nSymbols = 0; continue; }
          if(sym == '\n') continue;
          if(header == 1) continue;
          SEQ[nSymbols++] = sym;
          }
      if(P->aa)
        vr[idx_reads-1] = Target(KM->M[idx_model], D, idx_reads-1,
                          SEQ, nSymbols, Pos);
      else
        vr[idx_reads-1] = TargetAA(KM->M[idx_model], D, idx_reads-1,
                          SEQ, nSymbols, Pos);

      fclose(F);
      fclose(Pos);
      Free(name);
      Free(name2);

      S->average[idx_model] = sum / (double) FA->nReads;
      if(P->verbose) fprintf(stderr, "[>] Done!                           "
      "                    \n"); // SPACES ARE VALID!

      double sum1 = 0;
      for(idx_reads = 0 ; idx_reads < FA->nReads ; ++idx_reads)
        sum1 = sum1 + pow(((double) vr[idx_reads] - S->average[idx_model]), 2);
    
      S->variance[idx_model] = sum1 / (double) FA->nReads;
      S->stddev  [idx_model] = sqrt(S->variance[idx_model]);

      DeleteKModel(KM->M[idx_model]);
      }
  
    if(P->plots) PrintAveragePlot(S);
    }
  
  if(P->verbose) fprintf(stderr, "[>] Printing overall mRAWs statistics ...\n");
  PrintStats(S, P->min);
  if(P->verbose) fprintf(stderr, "[>] Done!\n"); 

  if(P->verbose) fprintf(stderr, "[>] Printing mRAWs element distribution ...\n");
  if(P->aa)
    PrintDist(D, P->min);
  else
    PrintDistAA(D, AL, P->min);
  if(P->verbose) fprintf(stderr, "[>] Done!\n"); 

  if(P->aa)
    {
    if(P->verbose)
      fprintf(stderr, "[>] Printing mRAWs nucleotide %% ...\n");
    PrintCG(D, P->min, P->plots, P->verbose);
    if(P->verbose) fprintf(stderr, "[>] Done!\n");
    }
  else
    {
    if(P->verbose)
      fprintf(stderr, "[>] Printing mRAWs aminoacids %% ...\n");
    PrintDistAAPerc(D, AL, P->min, P->plots, P->verbose);
    if(P->verbose) fprintf(stderr, "[>] Done!\n");
    }

/*
  if(P->aa){
    if(P->gc_prof){
      if(P->verbose)
        fprintf(stderr, "[>] Printing CG profiles ...\n");
      for(x = 0 ; x < FA->nReads ; ++x)
        PrintCGProfile(x);
      AverageCGPlot(FA->max_nsym);
      if(P->verbose) fprintf(stderr, "[>] Done!\n"); 
      }
    }

  if(P->plots)
    {
    if(P->verbose)
      fprintf(stderr, "[>] Printing SH plot ...\n");
    PrintGnuplot(P->min, P->min + P->nKmers - 1, FA->max_nsym);
    if(P->verbose)
      {
      fprintf(stderr, "[>] Shell plot in kplot.sh\n");
      fprintf(stderr, "[>] Done!\n");
      }
    }
*/

  RemoveStats(S);
  if(P->aa) RemoveDist(D);

  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

