#include <stdio.h>
#include <stdlib.h>

#include "msg.h"
#include "keys.h"
#include "colors.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintWarning(char *s)
  {
  fprintf(stderr, "%s[x] Error: %s%s\n", error_color, s, normal_color);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintMessage(char *s)
  {
  fprintf(stderr, "%s[>] %s%s\n", normal_color, s, normal_color);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintVersion(void)
  {
  fprintf(stdout, "v%u.%u\n", VERSION, RELEASE);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintMenu(void)
  {
  fprintf(stderr,
  "                                                        \n"
  "Program: %s [ Alignment-free and spatial-temporal analysis \n"
  "                  toolkit for large-scale multi-FASTA data ]\n" 
  "Version: %u.%u                                          \n"
  "                                                        \n"
  "Usage: %s <command> [options] < <file>                  \n"
  "                                                        \n"
  "Commands:                                              \n"
  "      %-10s   Moving average filter of a column float\n"
  "                   CSV file (the column to use is a parameter).\n"
  "      %-10s   Filters FASTA reads by characteristics:\n"
  "                   alphabet, completeness, length, CG quantity,\n"
  "                   multiple string patterns and pattern absence.\n"
  "      %-10s   Computes the alphabet frequencies for each \n"
  "                   FASTA read (it enables alphabet filtering).\n"
  "      %-10s   Computes Normalized Compression (NC) for all\n"
  "                   FASTA reads according to a compression level.\n"
  "      %-10s   Computes Normalized Compression Distance (NCD)\n"
  "                   for all FASTA reads according to a reference.\n"
  "      %-10s   Computes Relative Absent Words (RAWs) with\n"
  "                   CG quantity estimation for all RAWs.\n"
  "                                                      \n"
  "Help: %s <command> -h for accessing each command menu.\n"
  "                                                      \n",
  PNAME, VERSION, RELEASE, PNAME, LT_KEYS[1].key, LT_KEYS[2].key, 
  LT_KEYS[3].key, LT_KEYS[4].key, LT_KEYS[5].key, LT_KEYS[6].key, 
  PNAME);
  return;
  }
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Moving Average

void PrintMenuMA(void)
  {
  fprintf(stderr,
  "NAME                                                                    \n"
  "      %s %s                                                    \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      filters a float column file using a moving average.               \n"
  "                                                                        \n"
  "PARAMETERS                                                              \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu),                                     \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information),                             \n"
  "                                                                        \n"
  "      -i,  --ignore-first-line                                          \n"
  "           ignores first line,                                          \n"
  "                                                                        \n"
  "      -p,  --position                                                   \n"
  "           show index position for each entry,                          \n"
  "                                                                        \n"
  "      -w [INT],  --window [INT]                                         \n"
  "           window size (default: 47),                                   \n"
  "                                                                        \n"
  "      -c [INT],  --column [INT]                                         \n"
  "           column to filter (e.g. first column: -c 1),                  \n"
  "                                                                        \n"
  "      < [FILE]                                                          \n"
  "           Input column file target (e.g. data.txt) -- MANDATORY.       \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      %s %s [OPTION]... < [FILE]            \n"
  "                                                                        \n"
  "EXAMPLE                                                                 \n"
  "      %s %s -v -w %u < data.txt             \n"
  "                                                                        \n",
  PNAME, LT_KEYS[1].key, PNAME, LT_KEYS[1].key, PNAME, LT_KEYS[1].key,
  DEF_MA_WINDOW);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Filter by characteristics

void PrintMenuFC(void)
  {
  fprintf(stderr,
  "NAME                                                                    \n"
  "      %s %s                                                     \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      Filters a multi-FASTA file by specific characteristics:           \n"
  "      alphabet, completeness, length, CG quantity, string patterns.     \n"
  "                                                                        \n"
  "PARAMETERS                                                              \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu),                                     \n"
  "                                                                        \n"
  "      -V,  --version                                                    \n"
  "           display program and version information,                     \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information),                             \n"
  "                                                                        \n"
  "      -c,  --complete                                                   \n"
  "           considers only genomes with matching alphabet,               \n"
  "                                                                        \n"
  "      -a [STRING],  --alphabet [STRING]                                 \n"
  "           alphabet to consider (Default: ACGT),                        \n"
  "                                                                        \n"
  "      -min [INT],  --minimum [INT]                                      \n"
  "           minimum sequence size,                                       \n"
  "                                                                        \n"
  "      -max [INT],  --maximum [INT]                                      \n"
  "           maximum sequence size,                                       \n"
  "                                                                        \n"
  "      -ncg [FLOAT],  --cg-minimum [FLOAT]                               \n"
  "           minimum CG quantity (between 0.0 and 1.0),                   \n"
  "                                                                        \n"
  "      -mcg [FLOAT],  --cg-maximum [FLOAT]                               \n"
  "           maximum CG quantity (between 0.0 and 1.0),                   \n"
  "                                                                        \n"
  "      -p [STRING],  --pattern [STRING]                                  \n"
  "           considers reads with pattern(s) in the header,               \n"
  "           repeat flag for multiple patterns (case insensitive),        \n"
  "                                                                        \n"
  "      -i [STRING],  --ignore [STRING]                                   \n"
  "           ignores reads with pattern(s) in the header,                 \n"
  "           repeat flag for multiple ignores (case insensitive),         \n"
  "                                                                        \n"
  "      < [FILE]                                                          \n"
  "           Input FASTA target (e.g. SARS-CoV-2.fa) -- MANDATORY.        \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      %s %s [OPTION]... < [FILE]                              \n"
  "                                                                        \n"
  "EXAMPLE                                                                 \n"
  "      %s %s -a ACGT -min 29200 -p SARS < viruses.fa           \n"
  "                                                                        \n",
  PNAME, LT_KEYS[2].key, PNAME, LT_KEYS[2].key, PNAME, LT_KEYS[2].key);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Frequencies

void PrintMenuAF(void)
  {
  fprintf(stderr,
  "NAME                                                                    \n"
  "      %s %s                                                   \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      computes the alphabet frequencies for each FASTA read.            \n"
  "                                                                        \n"
  "PARAMETERS                                                              \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu),                                     \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information),                             \n"
  "                                                                        \n"
  "      -p,  --first-line                                                 \n"
  "           print first line with symbols,                               \n"
  "                                                                        \n"
  "      -a [STRING],  --alphabet [STRING]                                 \n"
  "           alphabet to consider (Default: ACGT),                        \n"
  "                                                                        \n"
  "      < [FILE]                                                          \n"
  "           Input FASTA target (e.g. SARS-CoV-2.fa) -- MANDATORY.        \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      %s %s [OPTION]... < [FILE]                              \n"
  "                                                                        \n"
  "EXAMPLE                                                                 \n"
  "      %s %s < viruses.fa                                      \n"
  "                                                                        \n",
  PNAME, LT_KEYS[3].key, PNAME, LT_KEYS[3].key, PNAME, LT_KEYS[3].key);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintModels(void){
  fprintf(stderr,
  "                                                               \n"
  "-m [C]:[D]:[R]:[I]:[H]:[G]/[S]:[E]:[A]                         \n"
  "                                                               \n"
  "Template of a target context model.                            \n"
  "                                                               \n"
  "Parameters:                                                    \n"
  "     [C]: (integer [1;20]) order size of the regular context   \n"
  "          model. Higher values use more RAM but, usually, are  \n"
  "          related to a better compression score.               \n"
  "     [D]: (integer [1;5000]) denominator to build alpha, which \n"
  "          is a parameter estimator. Alpha is given by 1/[D].   \n"
  "          Higher values are usually used with higher [C],      \n"
  "          and related to confiant bets. When [D] is one,       \n"
  "          the probabilities assume a Laplacian distribution.   \n"
  "     [R]: (integer [0;99999999] memory model. The 0 uses the   \n"
  "          full memory.                                         \n"
  "     [I]: (integer {0,1,2}) number to define if a sub-program  \n"
  "          which addresses the specific properties of DNA       \n"
  "          sequences (Inverted repeats) is used or not. The     \n"
  "          number 2 turns ON this sub-program without the       \n"
  "          regular context model (only inverted repeats). The   \n"
  "          number 1 turns ON the sub-program using at the same  \n"
  "          time the regular context model. The number 0 does    \n"
  "          not contemple its use (Inverted repeats OFF). The    \n"
  "          use of this sub-program increases the necessary time \n"
  "          to compress but it does not affect the RAM.          \n"
  "     [H]: (integer [1;254]) size of the cache-hash for deeper  \n"
  "          context models, namely for [C] > 14. When the        \n"
  "          [C] <= 14 use, for example, 1 as a default. The      \n"
  "          RAM is highly dependent of this value (higher value  \n"
  "          stand for higher RAM).                               \n"
  "     [G]: (real [0;1)) real number to define gamma. This value \n"
  "          represents the decayment forgetting factor of the    \n"
  "          regular context model in definition.                 \n"
  "     [S]: (integer [0;20]) maximum number of editions allowed  \n"
  "          to use a substitutional tolerant model with the same \n"
  "          memory model of the regular context model with       \n"
  "          order size equal to [C]. The value 0 stands for      \n"
  "          turning the tolerant context model off. When the     \n"
  "          model is on, it pauses when the number of editions   \n"
  "          is higher that [C], while it is turned on when       \n"
  "          a complete match of size [C] is seen again. This     \n"
  "          is probabilistic-algorithmic model very usefull to   \n"
  "          handle the high substitutional nature of genomic     \n"
  "          sequences. When [S] > 0, the compressor used more    \n"
  "          processing time, but uses the same RAM and, usually, \n"
  "          achieves a substantial higher compression ratio. The \n"
  "          impact of this model is usually only noticed for     \n"
  "          [C] >= 14.                                           \n"
  "     [E]: (integer [1;5000]) denominator to build alpha for    \n"
  "          substitutional tolerant context model. It is         \n"
  "          analogous to [D], however to be only used in the     \n"
  "          probabilistic model for computing the statistics of  \n"
  "          the substitutional tolerant context model.           \n"
  "     [A]: (real [0;1)) real number to define gamma. This value \n"
  "          represents the decayment forgetting factor of the    \n"
  "          substitutional tolerant context model in definition. \n"
  "          Its definition and use is analogus to [G].           \n"
  "                                                               \n");
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Normalized Compression

void PrintMenuNC(void){
  fprintf(stderr,
  "NAME                                                                \n"
  "      %s %s                                             \n"
  "                                                                    \n"
  "DESCRIPTION                                                         \n"
  "      Normalized compression estimation of multiple FASTA reads.    \n"
  "                                                                    \n"
  "PARAMETERS                                                          \n"
  "                                                                    \n"
  "      -h,  --help                                                   \n"
  "           usage guide (help menu),                                 \n"
  "                                                                    \n"
  "      -v,  --verbose                                                \n"
  "           verbose mode (more information),                         \n"
  "                                                                    \n"
  "      -d,  --dna                                                    \n"
  "           considers exclusively DNA alphabet {A,C,G,T},            \n"
  "           it also provides inverted repeats models,                \n"
  "           flag absence considers inversions (without complements), \n"
  "                                                                    \n"
  "      -p,  --show-parameters                                        \n"
  "           show parameters of the models for optimization,          \n"
  "                                                                    \n"
  "      -s,  --show-levels                                            \n"
  "           show pre-computed compression levels (parameters),       \n"
  "                                                                    \n"
  "      -t [INT],  --threads [INT]                                    \n"
  "           maximum number of threads to compute it,                 \n"
  "                                                                    \n"
  "      -l [INT],  --level [INT]                                      \n"
  "           compression level (integer),                             \n"
  "           it defines compressibility in balance with computational \n"
  "           resources (RAM & time), use -s for levels perception,    \n"
  "                                                                    \n"
  "      [FILE]                                                        \n"
  "           input sequence filename (to analyze) -- MANDATORY,       \n"
  "           FASTA file for the analysis (last argument).             \n"
  "                                                                    \n"
  "SYNOPSIS                                                            \n"
  "      %s %s [OPTION]... [FILE]                           \n"
  "                                                                    \n"
  "EXAMPLE                                                             \n"
  "      %s %s -v -t 8 -m 11:50:0:1:0:0.9/0:0:0 seqs.fa     \n"
  "                                                                    \n",
  PNAME, LT_KEYS[4].key, PNAME, LT_KEYS[4].key, PNAME, LT_KEYS[4].key);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Normalized Compression Distance

void PrintMenuNCD(void){
  fprintf(stderr,
  "NAME                                                                \n"
  "      %s %s                                             \n"
  "                                                                    \n"
  "DESCRIPTION                                                         \n"
  "      Normalized Compression Distance (NCD) between a given         \n" 
  "      reference sequence and each read in a multi-FASTA file.       \n"
  "                                                                    \n"
  "PARAMETERS                                                          \n"
  "                                                                    \n"
  "      -h,  --help                                                   \n"
  "           usage guide (help menu),                                 \n"
  "                                                                    \n"
  "      -v,  --verbose                                                \n"
  "           verbose mode (more information),                         \n"
  "                                                                    \n"
  "      -d,  --dna                                                    \n"
  "           considers exclusively DNA alphabet {A,C,G,T},            \n"
  "           it also provides inverted repeats models,                \n"
  "           flag absence considers inversions (without complements), \n"
  "                                                                    \n"
  "      -p,  --show-parameters                                        \n"
  "           show parameters of the models for optimization,          \n"
  "                                                                    \n"
  "      -s,  --show-levels                                            \n"
  "           show pre-computed compression levels (parameters),       \n"
  "                                                                    \n"
  "      -t [INT],  --threads [INT]                                    \n"
  "           maximum number of threads to compute it,                 \n"
  "                                                                    \n"
  "      -l [INT],  --level [INT]                                      \n"
  "           compression level (integer),                             \n"
  "           it defines compressibility in balance with computational \n"
  "           resources (RAM & time), use -s for levels perception,    \n"
  "                                                                    \n"
  "      -r [FILE],  --reference [FILE]                                \n"
  "           reference sequence to calculate the NCD to the reads,    \n"
  "                                                                    \n"
  "      [FILE]                                                        \n"
  "           input sequence filename (to analyze) -- MANDATORY,       \n"
  "           FASTA file for the analysis (last argument).             \n"
  "                                                                    \n"
  "SYNOPSIS                                                            \n"
  "      %s %s [OPTION]... [FILE]                           \n"
  "                                                                    \n"
  "EXAMPLE                                                             \n"
  "      %s %s -v -m 11:50:0:1:0:0.9/0:0:0 -r ref.fa seqs.fa \n"
  "                                                                    \n",
  PNAME, LT_KEYS[5].key, PNAME, LT_KEYS[5].key, PNAME, LT_KEYS[5].key);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RAWs

void PrintMenuRW(void){
  fprintf(stderr,
  "NAME                                                                \n"
  "      %s %s                                   \n"
  "                                                                    \n"
  "DESCRIPTION                                                         \n"
  "      Computation of minimal Relative Absent Words (mRAWs).         \n"
  "                                                                    \n"
  "PARAMETERS                                                          \n"
  "                                                                    \n"
  "      -h,  --help                                                   \n"
  "           usage guide (help menu),                                 \n"
  "                                                                    \n"
  "      -f,  --force                                                  \n"
  "           force mode (overwrites old files)                        \n"
  "                                                                    \n"
  "      -v,  --verbose                                                \n"
  "           verbose mode (more information)                          \n"
  "                                                                    \n"
  "      -vv, --very-verbose                                           \n"
  "           very verbose mode (much more information)                \n"
  "                                                                    \n"
  "      -a,  --aminoacids                                             \n"
  "           use amino acids/proteins models                          \n"
  "                                                                    \n"
  "      -t,  --threads                                                \n"
  "           does NOT use threads if flag is set (slower)             \n"
  "                                                                    \n"
  "      -i,  --ignore-ir                                              \n"
  "           does NOT use inverted repeats if flag is set             \n"
  "                                                                    \n"
  "      -o,  --stdout                                                 \n"
  "           write overall statistics to standard output              \n"
  "                                                                    \n"
  "      -p,  --plots                                                  \n"
  "           print Shell code to generate plots (gnuplot)             \n"
  "                                                                    \n"
  "      -min [NUMBER],  --minimum [NUMBER]                            \n"
  "           k-mer minimum size (usually 10)                          \n"
  "                                                                    \n"
  "      -max [NUMBER],  --maximum [NUMBER]                            \n"
  "           k-mer maximum size (usually 16)                          \n"
  "                                                                    \n"
  "      [FILE]                                                        \n"
  "           Input host FASTA (e.g. human) -- MANDATORY.              \n"
  "           This content will be loaded in the models.               \n"
  "                                                                    \n"
  "      [FILE]                                                        \n"
  "           Input parasite FASTA  (e.g. SARS-CoV-2) -- MANDATORY.    \n"
  "           The mRAWs will be mapped on this content file.           \n"
  "                                                                    \n"
  "SYNOPSIS                                                            \n"
  "      %s %s [OPTION]... > output.fa      \n"
  "                                                                    \n"
  "EXAMPLE                                                             \n"
  "      %s %s -v -min 11 -max 16 human.fa SARS-CoV2.fa\n"
  "                                                                    \n",
  PNAME, LT_KEYS[6].key, PNAME, LT_KEYS[6].key, PNAME, LT_KEYS[6].key);
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

