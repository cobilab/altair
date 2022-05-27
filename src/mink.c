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
  FILE *IN;
  FILE *OUT = NULL;
  char tmp[65535];
  unsigned x;

  if(argc < 4) {
    printf("Usage: mink <first> <second> <FILE>\n");
    return 1;
    }

  printf("Starting mink ...\n");
  IN  = Fopen(argv[argc-1], "r");
  OUT = Fopen(Cat(argv[argc-1], ".mink"), "w");
    
  while(fscanf(IN, "%u\t%s\n", &x, tmp) == 2)
      fprintf(OUT, "%u\t%u\t%u\n", atoi(argv[argc-3]), x, atoi(argv[argc-2])); 

  fclose(IN);
  fclose(OUT);
  printf("Done all!\n");
  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
