#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *Cat(char *a, char *b){
  char *base = (char *) malloc(strlen(a) + strlen(b) + 1);
  strcpy(base, a);
  strcat(base, b);
  return base;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE *Fopen(const char *p, const char *m){
  FILE *F = fopen(p, m);
  if(F == NULL){
    fprintf(stderr, "Error opening: %s (mode %s). Does it exist?\n", p, m);
    exit(1);
    }
  return F;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main(int argc, char *argv[]){
  FILE *IN  = NULL;
  FILE *OUT = NULL;
  unsigned x, y, strain, pos, order, **matrix, nSeqs, length, maxk;
  
  if(argc < 4) {
    printf("Usage: rebat <nSeqs> <length> <maxk> <FILE>\n");
    return 1;
    }

  printf("Starting rebat ...\n");
  IN     = Fopen(argv[argc-1], "r");
  OUT    = Fopen("out", "w");
  length = atoi(argv[argc-3]);
  maxk   = atoi(argv[argc-2]) + 1;
  nSeqs  = atoi(argv[argc-4]) + 1;

  matrix = (unsigned **) malloc(nSeqs * sizeof(unsigned *));
  for(x = 1 ; x < nSeqs ; ++x){
    matrix[x] = (unsigned *) malloc(length * sizeof(unsigned));
    for(y = 1 ; y < length ; ++y)
      matrix[x][y] = maxk;
    }

  while(fscanf(IN, "%u\t%u\t%u\n", &strain, &pos, &order) == 3){
    if(matrix[strain][pos] > order)
      matrix[strain][pos] = order;  
    }

  fprintf(OUT, "# X	Y	Z\n");
  for(x = 1 ; x < nSeqs ; ++x)
    for(y = 1 ; y < length ; ++y)
      if(matrix[x][y] != maxk)
        fprintf(OUT, "%u\t%u\t%u\n", y, x, matrix[x][y]);

  printf("Done!\n");
  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
