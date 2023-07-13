#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>

#include "ncd.h"
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

NCD_PARAMETERS *NP;        // FOR THREAD SHARING
FASTA_READS    *NFA;       // SHARE PARAMETERS FOR THREADS
CMODEL         **cModels;  // SHARED FOR THREADING

// - - - - - - - - - - C O M P R E S S O R   D N A - - - - - - - - - - - - - -

uint64_t CompressTargetReadAndLoadModels(FILE *F)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CACHE       *C;
  double      bits = 0, bps = 0;

  // FORCE DNA
  NP->nSym = 4;

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM            = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
        PM[n]   = CreatePModel      (NP->nSym);
  MX            = CreatePModel      (NP->nSym);
  PT            = CreateFloatPModel (NP->nSym);
  WM            = CreateWeightModel (totModels);

  char *readBUF;
  readBUF  = (uint8_t *) Calloc(DEF_BUF_SIZE, sizeof(uint8_t));
  symbBUF  = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(NP->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  compressed = 0;
  uint8_t header = 1;
  while((k = fread(readBUF, 1, DEF_BUF_SIZE, F)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = readBUF[idxPos];

      if(sym == '>'){ header = 1; continue; }
      if(sym == '\n' && header == 1){ header = 0; continue; }
      if(sym == '\n') continue;
      if(header == 1) continue;

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
        {
        bits += 2;
        ++compressed;
        continue;
        }

      symbBUF->buf[symbBUF->idx] = sym = DNASymToNum(sym);
      memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

      n = 0;
      pos = &symbBUF->buf[symbBUF->idx-1];
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
        {
        CMODEL *CM = cModels[cModel];
        if(CM->memory != 0)
          {
          if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
          else                            CM->M.pos++;
          }
        }

      RenormalizeWeights(WM);

      for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
        if(cModels[cModel]->edits != 0)
          UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

      UpdateCBuffer(symbBUF);
      ++compressed;
      }

  Free(MX);
  
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(PM[n]);  
  Free(PM);
  Free(PT);
  Free(readBUF);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }

// - - - - - - - - - - C O M P R E S S O R   D N A - - - - - - - - - - - - - -

uint64_t CompressTargetReadAndLoadModelsAA(FILE *F)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CACHE       *C;
  double      bits = 0, bps = 0;

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM            = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
        PM[n]   = CreatePModel      (NP->nSym);
  MX            = CreatePModel      (NP->nSym);
  PT            = CreateFloatPModel (NP->nSym);
  WM            = CreateWeightModel (totModels);

  char *readBUF;
  readBUF  = (uint8_t *) Calloc(DEF_BUF_SIZE, sizeof(uint8_t));
  symbBUF  = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(NP->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  uint8_t header = 1;
  compressed = 0;
  while((k = fread(readBUF, 1, DEF_BUF_SIZE, F)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = readBUF[idxPos];

      if(sym == '>'){ header = 1; continue; }
      if(sym == '\n' && header == 1){ header = 0; continue; }
      if(sym == '\n') continue;
      if(header == 1) continue;

      symbBUF->buf[symbBUF->idx] = sym = NP->A->revMap[sym];
      memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

      n = 0;
      pos = &symbBUF->buf[symbBUF->idx-1];
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

      for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
        {
        CMODEL *CM = cModels[cModel];
        if(CM->memory != 0)
          {
          if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
          else                            CM->M.pos++;
          }
        }

      RenormalizeWeights(WM);

      for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
        if(cModels[cModel]->edits != 0)
          UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

      UpdateCBuffer(symbBUF);
      ++compressed;
      }

  Free(MX);

  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(PM[n]);
  Free(PM);
  Free(PT);
  Free(readBUF);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }



//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - C O M P R E S S O R   D N A - - - - - - - - - - - - - -

uint64_t CompressTargetRead_NCD(char *SEQ, uint64_t length)
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
  NP->nSym = 4;
  
  // EXTRA MODELS DERIVED FROM EDITS
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM            = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
        PM[n]   = CreatePModel      (NP->nSym);
  MX            = CreatePModel      (NP->nSym);
  PT            = CreateFloatPModel (NP->nSym);
  WM            = CreateWeightModel (totModels);

  symbBUF  = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  cModels = (CMODEL **) Malloc(NP->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < NP->nModels ; ++n)
    cModels[n] = CreateCModel(NP->model[n].ctx, NP->model[n].den,
    TARGET, NP->model[n].edits, NP->model[n].eDen, NP->nSym, 
    NP->model[n].gamma, NP->model[n].eGamma, NP->model[n].ir, 
    NP->model[n].eIr, NP->model[n].memory);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(NP->model[n].edits != 0)
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
    memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
      if(cModels[cModel]->edits != 0)
        UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);
  for(n = 0 ; n < NP->nModels ; ++n)
     RemoveCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    Free(PM[n]->freqs);
    Free(PM[n]);
    }
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - C O M P R E S S O R   A A - - - - - - - - - - - - - - -

uint64_t CompressTargetReadAA_NCD(char *SEQ, uint64_t length)
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
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM      = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
    PM[n] = CreatePModel      (NP->nSym);
  MX      = CreatePModel      (NP->nSym);
  PT      = CreateFloatPModel (NP->nSym);
  WM      = CreateWeightModel (totModels);

  symbBUF = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  cModels = (CMODEL **) Malloc(NP->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < NP->nModels ; ++n)
    cModels[n] = CreateCModel(NP->model[n].ctx, NP->model[n].den,
    TARGET, NP->model[n].edits, NP->model[n].eDen, NP->nSym, 
    NP->model[n].gamma, NP->model[n].eGamma, NP->model[n].ir, 
    NP->model[n].eIr, NP->model[n].memory);
  
  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(NP->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  compressed = 0;
  for(idxPos = 0 ; idxPos < length ; ++idxPos)
    {
    sym = SEQ[idxPos];
    symbBUF->buf[symbBUF->idx] = sym = NP->A->revMap[sym];
      
    memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    ComputeMXProbs(PT, MX, NP->nSym);
    bits += (bps = PModelNats(MX, sym) / M_LN2);
    CalcDecayment(WM, PM, sym);

    // ADD COUNTERS
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
      if(cModels[cModel]->edits != 0)
        UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);
  for(n = 0 ; n < NP->nModels ; ++n)
    RemoveCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(PM[n]);
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - C O M P R E S S O R   D N A  ->  C(X,Y):  T H E   Y   P A R T - - - -

uint64_t CompressTargetRead_Con(char *SEQ, uint64_t length)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CACHE       *C;
  double      bits = 0, bps = 0;

  // FORCE DNA
  NP->nSym = 4;

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM            = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
        PM[n]   = CreatePModel      (NP->nSym);
  MX            = CreatePModel      (NP->nSym);
  PT            = CreateFloatPModel (NP->nSym);
  WM            = CreateWeightModel (totModels);

  symbBUF  = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

  CMODEL **CM_clone;
  CM_clone = (CMODEL **) Malloc(NP->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < NP->nModels ; ++n)
    CM_clone[n] = CloneCModel(cModels[n]);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = CM_clone[n]->gamma;
    if(NP->model[n].edits != 0)
      WM->gamma[pIdx++] = CM_clone[n]->eGamma;
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
    memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = CM_clone[cModel];
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

    ComputeMXProbs(PT, MX, CM_clone[0]->nSym);

    bits += (bps = PModelNats(MX, sym) / M_LN2);
    CalcDecayment(WM, PM, sym);

    // ADD COUNTERS
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = CM_clone[cModel];
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
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = CM_clone[cModel];
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = CM_clone[cModel];
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = CM_clone[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
      if(CM_clone[cModel]->edits != 0)
        UpdateTolerantModel(CM_clone[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);

  for(n = 0 ; n < NP->nModels ; ++n)
     RemoveCModel(CM_clone[n]);
  Free(CM_clone);

  for(n = 0 ; n < totModels ; ++n){
    Free(PM[n]->freqs);
    Free(PM[n]);
    }
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - C O M P R E S S O R   A A  ->  C(X,Y):  T H E   Y   P A R T - - - -

uint64_t CompressTargetReadAA_Con(char *SEQ, uint64_t length)
  {
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     sym, irSym, *pos;
  PMODEL      **PM, *MX;
  FPMODEL     *PT;
  CMWEIGHT    *WM;
  CBUF        *symbBUF;
  CACHE       *C;
  double      bits = 0, bps = 0;

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = NP->nModels;
  for(n = 0 ; n < NP->nModels ; ++n)
    if(NP->model[n].edits != 0)
      totModels += 1;

  PM      = (PMODEL  **) Calloc(totModels, sizeof(PMODEL *));
  for(n = 0 ; n < totModels ; ++n)
    PM[n] = CreatePModel      (NP->nSym);
  MX      = CreatePModel      (NP->nSym);
  PT      = CreateFloatPModel (NP->nSym);
  WM      = CreateWeightModel (totModels);

  symbBUF = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);

/*
  CMODEL      **cModels;
  cModels = (CMODEL **) Malloc(NP->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < NP->nModels ; ++n)
    cModels[n] = CreateCModel(NP->model[n].ctx, NP->model[n].den,
    TARGET, NP->model[n].edits, NP->model[n].eDen, NP->nSym,
    NP->model[n].gamma, NP->model[n].eGamma, NP->model[n].ir,
    NP->model[n].eIr, NP->model[n].memory);
*/
// TODO: CLONE THE cMODELS

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < NP->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(NP->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->eGamma;
    }

  compressed = 0;
  for(idxPos = 0 ; idxPos < length ; ++idxPos)
    {
    sym = SEQ[idxPos];
    symbBUF->buf[symbBUF->idx] = sym = NP->A->revMap[sym];

    memset((void *)PT->freqs, 0, NP->nSym * sizeof(double));

    n = 0;
    pos = &symbBUF->buf[symbBUF->idx-1];
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    ComputeMXProbs(PT, MX, NP->nSym);
    bits += (bps = PModelNats(MX, sym) / M_LN2);
    CalcDecayment(WM, PM, sym);

    // ADD COUNTERS
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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
    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
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

    for(cModel = 0 ; cModel < NP->nModels ; ++cModel)
      {
      CMODEL *CM = cModels[cModel];
      if(CM->memory != 0)
        {
        if(CM->M.pos >= CM->M.size - 1) CM->M.pos = 0;
        else                            CM->M.pos++;
        }
      }

    RenormalizeWeights(WM);

    for(cModel = 0, n = 0 ; cModel < NP->nModels ; ++cModel, ++n)
      if(cModels[cModel]->edits != 0)
        UpdateTolerantModel(cModels[cModel]->TM, PM[++n], sym);

    UpdateCBuffer(symbBUF);
    ++compressed;
    }

  Free(MX);

// TODO: REMOVE THE CLONED MODELS
/*
  for(n = 0 ; n < NP->nModels ; ++n)
    RemoveCModel(cModels[n]);
*/

  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(PM[n]);
  Free(PM);
  Free(PT);
  RemoveCBuffer(symbBUF);

  return (uint64_t) bits;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -
/*
void *CompressThread(void *Thr)
  {
  THREADS *T = (THREADS *) Thr;
  if(NP->dna == 0) CompressTargetDNA(T[0]);
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
// - - - N O R M A L I Z E D   C O M P R E S S I O N   D I S T A N C E - - - -
//
void NormalizedCompressionDistance(NCD_PARAMETERS *MAP)
  {
  NP = *&MAP;

  CheckFileIsFASTA(NP->filename);
  if(NP->verbose) PrintMessage("Loading alphabet ...");
  NP->A = CreateAlphabet();
  LoadAlphabet(NP->A, NP->filename);
  NP->nSym = NP->A->cardinality;
  if(NP->verbose) PrintAlphabet(NP->A);
  if(NP->dna == 0)
    {
    NP->nSym = 4;
    if(NP->verbose)
      PrintMessage("Adapting alphabet to 4 symbols {ACGT} (flag --dna on)");
    }

  FILE *F = Fopen(NP->filename, "r");
  uint8_t buffer[BUFFER_SIZE], sym = 0, header = 1;
  uint32_t k, idx;
  uint64_t nSymbols = 0, idx_reads = 0;
 
  if(NP->verbose) fprintf(stderr, "[>] Analyzing FASTA reads...\n");
  NFA = CreateFastaReads();
  CountFastaReadsAndMax(NFA, NP->filename);
  if(NP->verbose)
    {
    fprintf(stderr, "[>] Number of FASTA reads: %"PRIu64"\n", NFA->nReads);
    fprintf(stderr, "[>] Maximum read length: %u\n", NFA->max_nsym);
    }

  uint64_t vr_conjoint [NFA->nReads+1];
  uint64_t vr_regular  [NFA->nReads+1];
  char *SEQ = (char *) Calloc(NFA->max_nsym + 1, sizeof(char));

  if(NP->verbose) 
    fprintf(stderr, "[>] Compressing %s reference ...\n", !NP->dna ? 
    "DNA" : "Aminoacids");

  uint32_t n;
  cModels = (CMODEL **) Malloc(NP->nModels * sizeof(CMODEL *));
  for(n = 0 ; n < NP->nModels ; ++n)
    cModels[n] = CreateCModel(NP->model[n].ctx, NP->model[n].den,
    TARGET, NP->model[n].edits, NP->model[n].eDen, NP->nSym,
    NP->model[n].gamma, NP->model[n].eGamma, NP->model[n].ir,
    NP->model[n].eIr, NP->model[n].memory);

  uint64_t ref_bits;
  FILE *RF = Fopen(NP->reference, "r");
  if(!NP->dna) ref_bits = CompressTargetReadAndLoadModels   (RF);
  else        ref_bits = CompressTargetReadAndLoadModelsAA (RF);
  fclose(RF);
    
  if(NP->verbose) fprintf(stderr, "[>] Reference: %"PRIu64" bits\n", ref_bits);
  
  if(NP->verbose) fprintf(stderr, "[>] Compressing %"PRIu64" %s reads ...\n", 
		 NFA->nReads, !NP->dna ? "DNA" : "Aminoacids");
  
  while((k = fread(buffer, 1, BUFFER_SIZE, F)))
    for(idx = 0 ; idx < k ; ++idx)
      {
      sym = buffer[idx];
      if(sym == '>')
        {
        header = 1;
        if(idx_reads > 0)
          {
          if(!NP->dna)
	    {
            vr_regular [idx_reads-1] = CompressTargetRead_NCD(SEQ, nSymbols);
            vr_conjoint[idx_reads-1] = CompressTargetRead_Con(SEQ, nSymbols);
	    }
          else
	    {
	    vr_regular [idx_reads-1] = CompressTargetReadAA_NCD(SEQ, nSymbols);
	    vr_conjoint[idx_reads-1] = CompressTargetReadAA_Con(SEQ, nSymbols);
	    }
        
	  CalcProgress(NFA->nReads, idx_reads);
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
 
  if(!NP->dna)
    {
    vr_regular [idx_reads-1] = CompressTargetRead_NCD(SEQ, nSymbols);
    vr_conjoint[idx_reads-1] = CompressTargetRead_Con(SEQ, nSymbols);
    }
  else
    {
    vr_regular [idx_reads-1] = CompressTargetReadAA_NCD(SEQ, nSymbols);
    vr_conjoint[idx_reads-1] = CompressTargetReadAA_Con(SEQ, nSymbols);
    }

  for(idx_reads = 0 ; idx_reads < NFA->nReads ; ++idx_reads)
    fprintf(stdout, "%"PRIu64"\t%lf\n", idx_reads+1, 
    ((double) (ref_bits + vr_conjoint[idx_reads]) - 
    min(vr_regular[idx_reads], ref_bits)) / 
    max(vr_regular[idx_reads], ref_bits));

  fclose(F);
  
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
