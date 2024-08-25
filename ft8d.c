#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>

#include "pffft.h"

typedef float complex complex_t;
typedef float real_t;

struct SYNC
{
  int i, j;
  real_t s;
  real_t *p;
};

typedef struct SYNC sync_t;

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define N 174
#define M 83

#define NSTP 128
#define NSPS 1600
#define NFFT 3200

#define NSSY 5
#define NFOS 5

#define NSYM (93 * NSSY)

#define NTOKENS 2063592
#define MAX22 4194304

#define MAXGRID 32400

uint8_t mn[N][3] =
{
  {15, 44, 72},
  {24, 50, 61},
  {32, 57, 77},
  {0, 43, 44},
  {1, 6, 60},
  {2, 5, 53},
  {3, 34, 47},
  {4, 12, 20},
  {7, 55, 78},
  {8, 63, 68},
  {9, 18, 65},
  {10, 35, 59},
  {11, 36, 57},
  {13, 31, 42},
  {14, 62, 79},
  {16, 27, 76},
  {17, 73, 82},
  {21, 52, 80},
  {22, 29, 33},
  {23, 30, 39},
  {25, 40, 75},
  {26, 56, 69},
  {28, 48, 64},
  {2, 37, 77},
  {4, 38, 81},
  {45, 49, 72},
  {50, 51, 73},
  {54, 70, 71},
  {43, 66, 71},
  {42, 67, 77},
  {0, 31, 58},
  {1, 5, 70},
  {3, 15, 53},
  {6, 64, 66},
  {7, 29, 41},
  {8, 21, 30},
  {9, 17, 75},
  {10, 22, 81},
  {11, 27, 60},
  {12, 51, 78},
  {13, 49, 50},
  {14, 80, 82},
  {16, 28, 59},
  {18, 32, 63},
  {19, 25, 72},
  {20, 33, 39},
  {23, 26, 76},
  {24, 54, 57},
  {34, 52, 65},
  {35, 47, 67},
  {36, 45, 74},
  {37, 44, 46},
  {38, 56, 68},
  {40, 55, 61},
  {19, 48, 52},
  {45, 51, 62},
  {44, 69, 74},
  {26, 34, 79},
  {0, 14, 29},
  {1, 67, 79},
  {2, 35, 50},
  {3, 27, 50},
  {4, 30, 55},
  {5, 19, 36},
  {6, 39, 81},
  {7, 59, 68},
  {8, 9, 48},
  {10, 43, 56},
  {11, 38, 58},
  {12, 23, 54},
  {13, 20, 64},
  {15, 70, 77},
  {16, 29, 75},
  {17, 24, 79},
  {18, 60, 82},
  {21, 37, 76},
  {22, 40, 49},
  {6, 25, 57},
  {28, 31, 80},
  {32, 39, 72},
  {17, 33, 47},
  {12, 41, 63},
  {4, 25, 42},
  {46, 68, 71},
  {53, 54, 69},
  {44, 61, 67},
  {9, 62, 66},
  {13, 65, 71},
  {21, 59, 73},
  {34, 38, 78},
  {0, 45, 63},
  {0, 23, 65},
  {1, 4, 69},
  {2, 30, 64},
  {3, 48, 57},
  {0, 3, 4},
  {5, 59, 66},
  {6, 31, 74},
  {7, 47, 81},
  {8, 34, 40},
  {9, 38, 61},
  {10, 13, 60},
  {11, 70, 73},
  {12, 22, 77},
  {10, 34, 54},
  {14, 15, 78},
  {6, 8, 15},
  {16, 53, 62},
  {17, 49, 56},
  {18, 29, 46},
  {19, 63, 79},
  {20, 27, 68},
  {21, 24, 42},
  {12, 21, 36},
  {1, 46, 50},
  {22, 53, 73},
  {25, 33, 71},
  {26, 35, 36},
  {20, 35, 62},
  {28, 39, 43},
  {18, 25, 56},
  {2, 45, 81},
  {13, 14, 57},
  {32, 51, 52},
  {29, 42, 51},
  {5, 8, 51},
  {26, 32, 64},
  {24, 68, 72},
  {37, 54, 82},
  {19, 38, 76},
  {17, 28, 55},
  {31, 47, 70},
  {41, 50, 58},
  {27, 43, 78},
  {33, 59, 61},
  {30, 44, 60},
  {45, 67, 76},
  {5, 23, 75},
  {7, 9, 77},
  {39, 40, 69},
  {16, 49, 52},
  {41, 65, 67},
  {3, 21, 71},
  {35, 63, 80},
  {12, 28, 46},
  {1, 7, 80},
  {55, 66, 72},
  {4, 37, 49},
  {11, 37, 63},
  {58, 71, 79},
  {2, 25, 78},
  {44, 75, 80},
  {0, 64, 73},
  {6, 17, 76},
  {10, 55, 58},
  {13, 38, 53},
  {15, 36, 65},
  {9, 27, 54},
  {14, 59, 69},
  {16, 24, 81},
  {19, 29, 30},
  {11, 66, 67},
  {22, 74, 79},
  {26, 31, 61},
  {23, 68, 74},
  {18, 20, 70},
  {33, 52, 60},
  {34, 45, 46},
  {32, 58, 75},
  {39, 42, 82},
  {40, 41, 62},
  {48, 74, 82},
  {19, 43, 47},
  {41, 48, 56}
};

int16_t nm[M][7] =
{
  {3, 30, 58, 90, 91, 95, 152},
  {4, 31, 59, 92, 114, 145, -1},
  {5, 23, 60, 93, 121, 150, -1},
  {6, 32, 61, 94, 95, 142, -1},
  {7, 24, 62, 82, 92, 95, 147},
  {5, 31, 63, 96, 125, 137, -1},
  {4, 33, 64, 77, 97, 106, 153},
  {8, 34, 65, 98, 138, 145, -1},
  {9, 35, 66, 99, 106, 125, -1},
  {10, 36, 66, 86, 100, 138, 157},
  {11, 37, 67, 101, 104, 154, -1},
  {12, 38, 68, 102, 148, 161, -1},
  {7, 39, 69, 81, 103, 113, 144},
  {13, 40, 70, 87, 101, 122, 155},
  {14, 41, 58, 105, 122, 158, -1},
  {0, 32, 71, 105, 106, 156, -1},
  {15, 42, 72, 107, 140, 159, -1},
  {16, 36, 73, 80, 108, 130, 153},
  {10, 43, 74, 109, 120, 165, -1},
  {44, 54, 63, 110, 129, 160, 172},
  {7, 45, 70, 111, 118, 165, -1},
  {17, 35, 75, 88, 112, 113, 142},
  {18, 37, 76, 103, 115, 162, -1},
  {19, 46, 69, 91, 137, 164, -1},
  {1, 47, 73, 112, 127, 159, -1},
  {20, 44, 77, 82, 116, 120, 150},
  {21, 46, 57, 117, 126, 163, -1},
  {15, 38, 61, 111, 133, 157, -1},
  {22, 42, 78, 119, 130, 144, -1},
  {18, 34, 58, 72, 109, 124, 160},
  {19, 35, 62, 93, 135, 160, -1},
  {13, 30, 78, 97, 131, 163, -1},
  {2, 43, 79, 123, 126, 168, -1},
  {18, 45, 80, 116, 134, 166, -1},
  {6, 48, 57, 89, 99, 104, 167},
  {11, 49, 60, 117, 118, 143, -1},
  {12, 50, 63, 113, 117, 156, -1},
  {23, 51, 75, 128, 147, 148, -1},
  {24, 52, 68, 89, 100, 129, 155},
  {19, 45, 64, 79, 119, 139, 169},
  {20, 53, 76, 99, 139, 170, -1},
  {34, 81, 132, 141, 170, 173, -1},
  {13, 29, 82, 112, 124, 169, -1},
  {3, 28, 67, 119, 133, 172, -1},
  {0, 3, 51, 56, 85, 135, 151},
  {25, 50, 55, 90, 121, 136, 167},
  {51, 83, 109, 114, 144, 167, -1},
  {6, 49, 80, 98, 131, 172, -1},
  {22, 54, 66, 94, 171, 173, -1},
  {25, 40, 76, 108, 140, 147, -1},
  {1, 26, 40, 60, 61, 114, 132},
  {26, 39, 55, 123, 124, 125, -1},
  {17, 48, 54, 123, 140, 166, -1},
  {5, 32, 84, 107, 115, 155, -1},
  {27, 47, 69, 84, 104, 128, 157},
  {8, 53, 62, 130, 146, 154, -1},
  {21, 52, 67, 108, 120, 173, -1},
  {2, 12, 47, 77, 94, 122, -1},
  {30, 68, 132, 149, 154, 168, -1},
  {11, 42, 65, 88, 96, 134, 158},
  {4, 38, 74, 101, 135, 166, -1},
  {1, 53, 85, 100, 134, 163, -1},
  {14, 55, 86, 107, 118, 170, -1},
  {9, 43, 81, 90, 110, 143, 148},
  {22, 33, 70, 93, 126, 152, -1},
  {10, 48, 87, 91, 141, 156, -1},
  {28, 33, 86, 96, 146, 161, -1},
  {29, 49, 59, 85, 136, 141, 161},
  {9, 52, 65, 83, 111, 127, 164},
  {21, 56, 84, 92, 139, 158, -1},
  {27, 31, 71, 102, 131, 165, -1},
  {27, 28, 83, 87, 116, 142, 149},
  {0, 25, 44, 79, 127, 146, -1},
  {16, 26, 88, 102, 115, 152, -1},
  {50, 56, 97, 162, 164, 171, -1},
  {20, 36, 72, 137, 151, 168, -1},
  {15, 46, 75, 129, 136, 153, -1},
  {2, 23, 29, 71, 103, 138, -1},
  {8, 39, 89, 105, 133, 150, -1},
  {14, 57, 59, 73, 110, 149, 162},
  {17, 41, 78, 143, 145, 151, -1},
  {24, 37, 64, 98, 121, 159, -1},
  {16, 41, 74, 128, 169, 171, -1}
};

uint8_t nrw[M] =
{
  7,6,6,6,7,6,7,6,6,7,6,6,7,7,6,6,
  6,7,6,7,6,7,6,6,6,7,6,6,6,7,6,6,
  6,6,7,6,6,6,7,7,6,6,6,6,7,7,6,6,
  6,6,7,6,6,6,7,6,6,6,6,7,6,6,6,7,
  6,6,6,7,7,6,6,7,6,6,6,6,6,6,6,7,
  6,6,6
};

int ncw = 3;

int costas[7] = {3, 1, 4, 0, 6, 5, 2};

int graymap[8] = {0, 1, 3, 2, 5, 6, 4, 7};

char c0[38] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/";
char c1[37] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
char c2[36] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
char c3[10] = "0123456789";
char c4[27] = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";

complex_t *signal;
real_t window[NSPS], *map, llr[N];
sync_t *list;
uint8_t message[N];

complex_t *in, *out;
PFFFT_Setup *setup;

void sync()
{
  int i, j, k, m, n, jmax;
  real_t sum[2], s, smax;

  for(i = 0; i < NSYM; ++i)
  {
    for(j = 0; j < NSPS; ++j)
    {
      in[j] = window[j] * signal[i * NSTP + j];
    }

    pffft_transform_ordered(setup, in, out, NULL, PFFFT_FORWARD);

    for(j = 0; j < NFFT; ++j)
    {
      k = j < NFFT / 2 ? j + NFFT / 2 : j - NFFT / 2;
      map[i * NFFT + k] = cabsf(out[j]);
    }
  }

  for(i = 0; i < NFFT; ++i)
  {
    jmax = 0;
    smax = 0;
    for(j = -7 * NSSY; j < 21 * NSSY; ++j)
    {
      memset(sum, 0, sizeof(sum));
      for(k = 0; k < 7; ++k)
      {
        m = j + (k + 36) * NSSY;
        for(n = 0; n < 8; ++n) sum[0] += map[m * NFFT + i + n * NFOS];
        sum[1] += map[m * NFFT + i + costas[k] * NFOS];
      }
      s = 7 * sum[1] / (sum[0] - sum[1]);
      if(smax < s)
      {
        jmax = j;
        smax = s;
      }
    }
    list[i].i = i;
    list[i].j = jmax;
    list[i].s = smax;
    list[i].p = map + jmax * NFFT + i;
  }

  for(i = 2; i < NFFT - 2; ++i)
  {
    if((list[i - 2].s > list[i].s && list[i - 1].s > list[i].s) || (list[i + 1].s > list[i].s && list[i + 2].s > list[i].s))
    {
      list[i].p = NULL;
    }
  }
}

real_t max(real_t a, real_t b, real_t c, real_t d)
{
  real_t x, y;
  x = a > b ? a : b;
  y = c > d ? c : d;
  return x > y ? x : y;
}

void process(sync_t *cand)
{
  int i, j, k, l;
  real_t s[8], sum[2], var, sig;

  for(i = 0; i < 2; ++i)
  {
    for(j = 0; j < 29; ++j)
    {
      k = (i * 36 + j + 7) * NSSY * NFFT;

      for(l = 0; l < 8; ++l)
      {
        s[l] = log10f(1e-32 + cand->p[k + graymap[l] * NFOS]);
      }

      k = i * 87 + j * 3;

      llr[k + 0] = max(s[4], s[5], s[6], s[7]) - max(s[0], s[1], s[2], s[3]);
      llr[k + 1] = max(s[2], s[3], s[6], s[7]) - max(s[0], s[1], s[4], s[5]);
      llr[k + 2] = max(s[1], s[3], s[5], s[7]) - max(s[0], s[2], s[4], s[6]);
    }
  }

  memset(sum, 0, sizeof(sum));
  for(i = 0; i < N; ++i)
  {
    sum[0] += llr[i];
    sum[1] += llr[i] * llr[i];
  }
  var = (sum[1] - sum[0] * sum[0] / N) / N;
  sig = var > 0 ? sqrtf(var) : sqrtf(sum[1] / N);
  sig /= 4;

  for(i = 0; i < N; ++i)
  {
    llr[i] /= sig;
  }
}

int check()
{
  int i, j;
  uint8_t poly[15] = {1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1 };
  uint8_t data[96] = {0};

  for(i = 0; i < 77; ++i) data[i] = message[i];

  for(i = 0; i < 82; ++i)
  {
    if(data[i])
    {
      for(j = 0; j < 15; ++j)
      {
        data[i + j] = (data[i + j] + poly[j]) % 2;
      }
    }
  }

  for(i = 0; i < 14; ++i)
  {
    if(data[82+i] != message[77 + i]) return 0;
  }

  return 1;
}

int decode(int iterations)
{
  int i, j, k, l, iter, ibj, ichk, current, previous, counter;
  real_t x, x2, tnm, tmn, tov[N][3], toc[M][7], zn[N];

  memset(tov, 0, sizeof(tov));

  for(i = 0; i < M; ++i)
  {
    for(j = 0; j < nrw[i]; ++j)
    {
      toc[i][j] = llr[nm[i][j]];
    }
  }

  counter = 0;
  previous = 0;
  for(iter = 0; iter < iterations; ++iter)
  {
    for(i = 0; i < N; ++i)
    {
      zn[i] = llr[i];
      for(j = 0; j < ncw; ++j) zn[i] += tov[i][j];
      message[i] = zn[i] > 0;
    }

    current = 0;
    for(i = 0; i < M; ++i)
    {
      l = 0;
      for(j = 0; j < nrw[i]; ++j) l += message[nm[i][j]];
      if(l % 2 > 0) ++current;
    }

    if(current == 0 && check()) return 1;

    if(iter > 0)
    {
      ++counter;
      if(current < previous) counter = 0;
      if(counter > 4 && iter > 9 && current > 15) return 0;
    }

    previous = current;

    for(i = 0; i < M; ++i)
    {
      for(j = 0; j < nrw[i]; ++j)
      {
        ibj = nm[i][j];
        tnm = llr[ibj];
        for(k = 0; k < ncw; ++k)
        {
          if(mn[ibj][k] != i)
          {
            tnm += tov[ibj][k];
          }
        }
        x = -tnm / 2;
        x2 = x * x;
        toc[i][j] = x > 3.64 ? 1 : x < -3.64 ? -1 : (945 + (105 + x2) * x2) * x / (945 + (420 + 15 * x2) * x2);
      }
    }

    for(i = 0; i < N; ++i)
    {
      for(j = 0; j < ncw; ++j)
      {
        ichk = mn[i][j];
        tmn = 1;
        for(k = 0; k < nrw[ichk]; ++k)
        {
          if(nm[ichk][k] != i)
          {
            tmn *= toc[ichk][k];
          }
        }
        x = -tmn;
        x2 = x * x;
        tov[i][j] = 2 * (945 - (735 - 64 * x2) * x2) * x / (945 - (1050 - 225 * x2) * x2);
      }
    }
  }

  return 0;
}

void trim(char *s)
{
  char *p = s;
  int l = strlen(p);

  while(isspace(p[l - 1])) p[--l] = 0;
  while(*p && isspace(*p)) ++p, --l;

  memmove(s, p, l + 1);
}

int unpack(char *call, char *grid)
{
  int i, n;
  uint64_t icall;
  uint16_t igrid;
  uint8_t i3, iflip;

  call[0] = 0;
  grid[0] = 0;

  i3 = 0;
  for(i = 0; i < 3; ++i) i3 |= message[76 - i] << i;

  if(i3 == 1 || i3 == 2)
  {
    icall = 0;
    for(i = 0; i < 28; ++i) icall |= message[56 - i] << i;

    igrid = 0;
    for(i = 0; i < 15; ++i) igrid |= message[73 - i] << i;

    if(igrid <= MAXGRID && igrid != 32373)
    {
      grid[4] = 0;
      n = igrid;
      grid[3] = '0' + (n % 10);
      n /= 10;
      grid[2] = '0' + (n % 10);
      n /= 10;
      grid[1] = 'A' + (n % 18);
      n /= 18;
      grid[0] = 'A' + (n % 18);
    }

    n = icall - NTOKENS - MAX22;

    if(n < 0) return 0;

    call[6] = 0;
    for(i = 0; i < 3; ++i)
    {
      call[5 - i] = c4[n % 27];
      n /= 27;
    }
    call[2] = c3[n % 10];
    n /= 10;
    call[1] = c2[n % 36];
    n /= 36;
    call[0] = c1[n % 37];
    trim(call);

    return 1;
  }

  if(i3 == 4)
  {
    iflip = message[70];

    if(iflip) return 0;

    icall = 0;
    for(i = 0; i < 58; ++i) icall |= message[69 - i] << i;
    call[11] = 0;
    for(i = 0; i < 11; ++i)
    {
      call[10 - i] = c0[icall % 38];
      icall /= 38;
    }
    trim(call);

    return 1;
  }

  return 0;
}

int snr(sync_t *cand)
{
  return floor(20.0 * log10f(1e-32 + cand->s) - 26 + 0.5);
}

int main(int argc, char **argv)
{
  FILE *fp;
  double dialfreq;
  int i, j, freq;
  sync_t *curr, *next, temp;
  char call[12], grid[5];
  real_t a[4] = {0.35875, 0.48829, 0.14128, 0.01168};

  if(argc != 2)
  {
    return EXIT_FAILURE;
  }

  if((fp = fopen(argv[1], "rb")) == NULL)
  {
    fprintf(stderr, "Cannot open input file %s.\n", argv[1]);
    return EXIT_FAILURE;
  }

  signal = malloc(sizeof(complex_t) * 60000);
  map = malloc(sizeof(real_t) * NSYM * NFFT);
  list = malloc(sizeof(sync_t) * NFFT);

  in = pffft_aligned_malloc(sizeof(complex_t) * NFFT);
  out = pffft_aligned_malloc(sizeof(complex_t) * NFFT);
  setup = pffft_new_setup(NFFT, PFFFT_COMPLEX);

  for(i = 0; i < NSPS; ++i)
  {
    window[i] = a[0] -
      a[1] * cosf(2.0 * M_PI * i / (NSPS - 1)) +
      a[2] * cosf(4.0 * M_PI * i / (NSPS - 1)) -
      a[3] * cosf(6.0 * M_PI * i / (NSPS - 1));
  }

  memset(in, 0, sizeof(complex_t) * NFFT);

  fread(&dialfreq, 1, 8, fp);

  for(i = 0; i < 4; ++i)
  {
    fread(signal, 1, 480000, fp);

    sync();

    for(j = 2; j < NFFT - 2; ++j)
    {
      curr = &list[j];
      next = &list[j + 1];

      if(curr->p == NULL || curr->s < 3.0) continue;

      if(next->p != NULL && next->s > curr->s)
      {
        temp = *curr;
        *curr = *next;
        *next = temp;
      }

      process(curr);

      if(!decode(30) || !unpack(call, grid)) continue;

      next->p = NULL;

      freq = floor(dialfreq + j * 4.0e3 / NFFT - 2000 + 0.5);
      printf("%1d %4d %4d %5.2f %3d %8d %6s %4s\n", i, curr->i, curr->j, curr->s, snr(curr), freq, call, grid);
    }
  }

  return EXIT_SUCCESS;
}
