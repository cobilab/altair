#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>

#include "nc.h"
#include "defs.h"
#include "threads.h"
#include "mem.h"
#include "msg.h"
#include "dna.h"
#include "buffer.h"
#include "levels.h"
#include "param.h"
#include "alphabet.h"
#include "common.h"
#include "cache.h"
#include "pmodels.h"
#include "strings.h"
#include "cm.h"

NC_PARAMETERS *P;  // FOR THREAD SHARING
FASTA_READS   *FA; // SHARE PARAMETERS FOR THREADS

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - C O M P R E S S O R   D N A - - - - - - - - - - - - - -

double CompressTargetRead(char *SEQ, uint64_t length)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CMODEL      **cModels;
  CACHE       *C;
  double      bits = 0, bps = 0;
 
  // FORCE DNA
  P->nSym = 4;
  
  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].edits != 0)
      totModels += 1;

  PM            = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
        PM[n]   = CreatePModel      (P->nSym);
  MX            = CreatePModel      (P->nSym);
  PT            = CreateFloatPModel (P->nSym);
  WM            = CreateWeightModel (totModels);

  symbBUF  = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  cModels = (CMODEL **) Malloc(P->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < P->nModels ; ++n)
    cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den,
    TARGET, P->model[n].edits, P->model[n].eDen, P->nSym, 
    P->model[n].gamma, P->model[n].eGamma, P->model[n].ir, 
    P->model[n].eIr, P->model[n].memory);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  compressed = 0;
  for(idxPos = 0 ; idxPos < length ; ++idxPos)
    {
    sym = SEQ[idxPos];

    // FINAL FILTERING DNA CONTENT
    if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
      {
      bits += 2;
      ++compressed;
      continue;
      }

    symbBUF->buf[symbBUF->idx] = sym = DNASymToNum(sym);
    memset((void *)PT->freqs, 0, P->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      GetPModelIdx(pos, CM);
      ComputePModel(CM, PM[n], CM->pModelIdx, CM->alphaDen);
      ComputeWeightedFreqs(WM->weight[n], PM[n], PT, CM->nSym);
      if(CM->edits != 0)
        {
        ++n;
        CM->TM->seq->buf[CM->TM->seq->idx] = sym;
        CM->TM->idx = GetPModelIdxCorr(CM->TM->seq->buf+
        CM->TM->seq->idx-1, CM, CM->TM->idx);
        ComputePModel(CM, PM[n], CM->TM->idx, CM->TM->den);
        ComputeWeightedFreqs(WM->weight[n], PM[n], PT, CM->nSym);
        }
      ++n;
      }

    ComputeMXProbs(PT, MX, cModels[0]->nSym);

    bits += (bps = PModelNats(MX, sym) / M_LN2);
    CalcDecayment(WM, PM, sym);

    // ADD COUNTERS
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      switch(CM->ir)
        {
        case 0:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        break;
        case 1:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, CM);
        UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
        break;
        case 2:
        irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, CM);
        UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
        break;
        default:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        break;
        }
      }

    // UPDATE INDEXES & SYM CACHE
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        uint32_t mem_pos = (CM->M.pos) % CM->M.size;
        switch(CM->ir)
          {
          case 0:
          CM->M.idx[mem_pos] = CM->pModelIdx;
          CM->M.sym[mem_pos] = sym;
          break;
          case 1:
          CM->M.idx[mem_pos]    = CM->pModelIdx;
          CM->M.sym[mem_pos]    = sym;
          CM->M.idx_ir[mem_pos] = CM->pModelIdxIR;
          CM->M.sym_ir[mem_pos] = irSym;
          break;
          case 2:
          CM->M.idx_ir[mem_pos] = CM->pModelIdxIR;
          CM->M.sym_ir[mem_pos] = irSym;
          break;
          default:
          PrintWarning("no store action!");
          exit(1);
          }
        }
      }

    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        uint32_t mps = (CM->M.pos + 1) % CM->M.size;
        if(compressed > CM->M.size)
          {
          switch(CM->ir)
            {
            case 0:
            UnUpdateCModelCounter(CM, CM->M.sym[mps], CM->M.idx[mps]);
            break;
            case 1:
            UnUpdateCModelCounter(CM, CM->M.sym   [mps], CM->M.idx   [mps]);
            UnUpdateCModelCounter(CM, CM->M.sym_ir[mps], CM->M.idx_ir[mps]);
            break;
            case 2:
            UnUpdateCModelCounter(CM, CM->M.sym_ir[mps], CM->M.idx_ir[mps]);
            break;
            default:
            PrintWarning("no update action!");
            exit(1);
            }
          }
        }
      }

    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < P->nModels ; ++cModel, ++n)
      if(cModels[cModel]->edits != 0)
        UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);
  for(n = 0 ; n < P->nModels ; ++n)
     RemoveCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    Free(PM[n]->freqs);
    Free(PM[n]);
    }
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return bits / (double) (compressed * log2(P->nSym));
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - C O M P R E S S O R   A A - - - - - - - - - - - - - - -
//

double CompressTargetReadAA(char *SEQ, uint64_t length)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CMODEL      **cModels;
  CACHE       *C;
  double      bits = 0, bps = 0;

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].edits != 0)
      totModels += 1;

  PM      = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
    PM[n] = CreatePModel      (P->nSym);
  MX      = CreatePModel      (P->nSym);
  PT      = CreateFloatPModel (P->nSym);
  WM      = CreateWeightModel (totModels);

  symbBUF = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  cModels = (CMODEL **) Malloc(P->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < P->nModels ; ++n)
    cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den,
    TARGET, P->model[n].edits, P->model[n].eDen, P->nSym, 
    P->model[n].gamma, P->model[n].eGamma, P->model[n].ir, 
    P->model[n].eIr, P->model[n].memory);
  
  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  compressed = 0;
  for(idxPos = 0 ; idxPos < length ; ++idxPos)
    {
    sym = SEQ[idxPos];
    symbBUF->buf[symbBUF->idx] = sym = P->A->revMap[sym];
      
    memset((void *)PT->freqs, 0, P->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      GetPModelIdx(pos, CM);
      ComputePModel(CM, PM[n], CM->pModelIdx, CM->alphaDen);
      ComputeWeightedFreqs(WM->weight[n], PM[n], PT, CM->nSym);
      if(CM->edits != 0)
        {
        ++n;
        CM->TM->seq->buf[CM->TM->seq->idx] = sym;
        CM->TM->idx = GetPModelIdxCorr(CM->TM->seq->buf+
        CM->TM->seq->idx-1, CM, CM->TM->idx);
        ComputePModel(CM, PM[n], CM->TM->idx, CM->TM->den);
        ComputeWeightedFreqs(WM->weight[n], PM[n], PT, CM->nSym);
        }
      ++n;
      }

    ComputeMXProbs(PT, MX, P->nSym);
    bits += (bps = PModelNats(MX, sym) / M_LN2);
    CalcDecayment(WM, PM, sym);

    // ADD COUNTERS
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      switch(CM->ir)
        {
        case 0:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        break;
        case 1:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        irSym = GetPModelIdxR(symbBUF->buf+symbBUF->idx, CM);
        UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
        break;
        case 2:
        irSym = GetPModelIdxR(symbBUF->buf+symbBUF->idx, CM);
        UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
        break;
        default:
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        break;
        }
      }

    // UPDATE INDEXES & SYM CACHE
    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        uint32_t mem_pos = (CM->M.pos) % CM->M.size;
        switch(CM->ir)
          {
          case 0:
          CM->M.idx[mem_pos] = CM->pModelIdx;
          CM->M.sym[mem_pos] = sym;
          break;
          case 1:
          CM->M.idx[mem_pos]    = CM->pModelIdx;
          CM->M.sym[mem_pos]    = sym;
          CM->M.idx_ir[mem_pos] = CM->pModelIdxIR;
          CM->M.sym_ir[mem_pos] = irSym;
          break;
          case 2:
          CM->M.idx_ir[mem_pos] = CM->pModelIdxIR;
          CM->M.sym_ir[mem_pos] = irSym;
          break;
          default:
          PrintWarning("no store action!");
          exit(1);
          }
        }
      }

    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        uint32_t mps = (CM->M.pos + 1) % CM->M.size;
        if(compressed > CM->M.size)
          {
          switch(CM->ir)
            {
            case 0:
            UnUpdateCModelCounter(CM, CM->M.sym[mps], CM->M.idx[mps]);
            break;
            case 1:
            UnUpdateCModelCounter(CM, CM->M.sym   [mps], CM->M.idx   [mps]);
            UnUpdateCModelCounter(CM, CM->M.sym_ir[mps], CM->M.idx_ir[mps]);
            break;
            case 2:
            UnUpdateCModelCounter(CM, CM->M.sym_ir[mps], CM->M.idx_ir[mps]);
            break;
            default:
            PrintWarning("no update action!");
            exit(1);
            }
          }
        }
      }

    for(cModel = 0 ; cModel < P->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < P->nModels ; ++cModel, ++n)
      if(cModels[cModel]->edits != 0)
        UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);
  for(n = 0 ; n < P->nModels ; ++n)
    RemoveCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(PM[n]);
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return bits / (double) (compressed * log2(P->nSym));
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -
/*
void *CompressThread(void *Thr)
  {
  THREADS *T = (THREADS *) Thr;
  if(P->dna == 0) CompressTargetDNA(T[0]);
  else            CompressTargetAA (T[0]);
  pthread_exit(NULL);
  }
*/
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - C O M P R E S S O R   M A I N - - - - - - - - - - - -
/*
void CompressAction(THREADS *T)
  {
  pthread_t t[2];
  uint32_t n;

  pthread_create(&(t[1]), NULL, CompressThread, (void *) &(T[0]));
  pthread_create(&(t[2]), NULL, CompressThread, (void *) &(T[1]));
    
  pthread_join(t[1], NULL);
  pthread_join(t[2], NULL);

  return;
  }
*/
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - N O R M A L I Z E D   C O M P R E S S I O N - - - - - - - -
//
void NormalizedCompression(NC_PARAMETERS *MAP)
  {
  P = *&MAP;

  CheckFileIsFASTA(P->filename);
  if(P->verbose) PrintMessage("Loading alphabet ...");
  P->A = CreateAlphabet();
  LoadAlphabet(P->A, P->filename);
  P->nSym = P->A->cardinality;
  if(P->verbose) PrintAlphabet(P->A);
  if(P->dna == 0)
    {
    P->nSym = 4;
    if(P->verbose)
      PrintMessage("Adapting alphabet to 4 symbols {ACGT} (flag --dna on)");
    }

  FILE *F = Fopen(P->filename, "r");
  uint8_t buffer[BUFFER_SIZE], sym = 0, header = 1;
  uint32_t k, idx;
  uint64_t nSymbols = 0, idx_reads = 0;
 
  if(P->verbose) fprintf(stderr, "[>] Analyzing FASTA reads...\n");
  FA = CreateFastaReads();
  CountFastaReadsAndMax(FA, P->filename);
  if(P->verbose)
    {
    fprintf(stderr, "[>] Number of FASTA reads: %"PRIu64"\n", FA->nReads);
    fprintf(stderr, "[>] Maximum read length: %u\n", FA->max_nsym);
    }
  if(P->verbose) fprintf(stderr, "[>] Done!\n");

  double vr[FA->nReads+1];
  char *SEQ = (char *) Calloc(FA->max_nsym + 1, sizeof(char));

  if(P->verbose) fprintf(stderr, "[>] Compressing %s ...\n", !P->dna ? 
		 "DNA" : "Aminoacids");
 
  char identifier_prefix[FA->nReads+1][HEADERS_PREFIX_SIZE+1];
  uint32_t idx_header = 0;

  while((k = fread(buffer, 1, BUFFER_SIZE, F)))
    for(idx = 0 ; idx < k ; ++idx)
      {
      sym = buffer[idx];
      if(sym == '>')
        {
        header = 1;
        if(idx_reads > 0)
          {
          if(!P->dna) vr[idx_reads-1] = CompressTargetRead   (SEQ, nSymbols);
          else        vr[idx_reads-1] = CompressTargetReadAA (SEQ, nSymbols);
          }
        ++idx_reads;
        continue;
        }
      if(sym == '\n' && header == 1)
        { 
	header = 0; 
	nSymbols = 0; 
	identifier_prefix[idx_reads-1][idx_header] = '\0';
	idx_header = 0;
       	continue; 
	}
      if(sym == '\n') continue;
      if(header == 1)
        {
	if(idx_header < HEADERS_PREFIX_SIZE)
	  identifier_prefix[idx_reads-1][idx_header++] = sym;
	continue;
        }
      SEQ[nSymbols++] = sym;
      }
 
  if(!P->dna) vr[idx_reads-1] = CompressTargetRead   (SEQ, nSymbols);
  else        vr[idx_reads-1] = CompressTargetReadAA (SEQ, nSymbols);
 
  for(idx_reads = 0 ; idx_reads < FA->nReads ; ++idx_reads)
    fprintf(stdout, "%"PRIu64"\t%lf\t%s\n", idx_reads+1, vr[idx_reads], 
    identifier_prefix[idx_reads]);

  fclose(F);
  
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
