#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "pffft.h"
#include "ldpc.h"

typedef float complex complex_t;
typedef float real_t;

#define NSTP 64
#define NSPS 1600
#define NFFT 3200

#define NSSY 10
#define NFOS 5

#define NSAMP 60000
#define NSYM (NSAMP / NSTP)

#define NTOKENS 2063592
#define MAX22 4194304

#define MAXGRID 32400

#define SYNC_MIN 2.0
#define SYNC_ITER_CAP 2.5

#define LDPC_ITER_LOW 20
#define LDPC_ITER_FULL 30

typedef struct
{
  int i, j, k;
  real_t r, s;
} sync_t;

static complex_t *signal;
static real_t window[NSPS], *map, llr[N];

static uint8_t nrw[M];
static int cpos[N][3];
static int p2[M][7];

static sync_t *list;
static uint8_t message[N];

static complex_t *buffer;
static PFFFT_Setup *setup;

static const int costas[7] = {3, 1, 4, 0, 6, 5, 2};
static const int graymap[8] = {0, 1, 3, 2, 5, 6, 4, 7};

static const char c0[38] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/";
static const char c1[37] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char c2[36] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char c3[10] = "0123456789";
static const char c4[27] = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/*
 * Sync stage: windowed FFTs and Costas search.
 */

/* Sum of the 8 oversampled bins of symbol m and the Costas tone at k. */
static inline void tone_sum(int i, int j, int k, real_t *sum, real_t *c)
{
  int m, n;
  real_t t;

  t = 0;
  m = j + (k + 36) * NSSY;

  for(n = 0; n < 8; ++n) t += map[m * NFFT + i + n * NFOS];

  *sum = t;
  *c = map[m * NFFT + i + costas[k] * NFOS];
}

static inline real_t sync_s(int i, int j)
{
  int k;
  real_t c, sum, s;

  s = 0;

  for(k = 0; k < 7; ++k)
  {
    tone_sum(i, j, k, &sum, &c);
    s += 8 * c / sum;
  }
  return s / 7;
}

static inline real_t sync_r(int i, int j)
{
  int k;
  real_t c, sum, r;

  r = 0;

  for(k = 0; k < 7; ++k)
  {
    tone_sum(i, j, k, &sum, &c);
    r += 7 * c / (sum - c);
  }
  return r / 7;
}

static void sync()
{
  int i, j, idx, len, jmax, jstp;
  real_t rmax, s, smax;
  real_t *mlo, *mhi;

  for(i = 0; i < NSYM; ++i)
  {
    memset(buffer + NSPS, 0, sizeof(complex_t) * (NFFT - NSPS));

    idx = i * NSTP;
    len = NSAMP - idx;
    if(len > NSPS) len = NSPS;
    for(j = 0; j < len; ++j)
    {
      buffer[j] = window[j] * signal[idx + j];
    }
    for(; j < NSPS; ++j)
    {
      buffer[j] = 0;
    }

    pffft_transform_ordered(setup, (float *)buffer, (float *)buffer, NULL, PFFFT_FORWARD);

    mlo = map + i * NFFT;
    mhi = mlo + NFFT / 2;
    for(j = 0; j < NFFT / 2; ++j)
    {
      mhi[j] = cabsf(buffer[j]);
      mlo[j] = cabsf(buffer[j + NFFT / 2]);
    }
  }

  /* Coarse-to-fine search of the start time for each frequency bin. */
  jstp = 5;

  for(i = 0; i < NFFT; ++i)
  {
    jmax = 0;
    smax = 0;

    for(j = -10 * NSSY; j < 25 * NSSY; j += jstp)
    {
      s = sync_s(i, j);

      if(s > smax)
      {
        jmax = j;
        smax = s;
      }
    }

    if(smax > 1.7)
    {
      for(j = jmax - jstp; j <= jmax + jstp; ++j)
      {
        s = sync_s(i, j);

        if(s > smax)
        {
          jmax = j;
          smax = s;
        }
      }
    }

    rmax = sync_r(i, jmax);

    list[i].i = i;
    list[i].j = jmax;
    list[i].k = 1;
    list[i].r = rmax;
    list[i].s = smax;
  }

  for(i = 2; i < NFFT - 2; ++i)
  {
    if((list[i - 2].s > list[i].s && list[i - 1].s > list[i].s) || (list[i + 1].s > list[i].s && list[i + 2].s > list[i].s))
    {
      list[i].k = 0;
    }
  }
}

static int count_costas(sync_t *cand)
{
  int i, j, m, n, nmax, result;
  real_t v, vmax;

  result = 0;

  for(i = 0; i < 3; ++i)
  {
    for(j = 0; j < 7; ++j)
    {
      m = cand->j + (i * 36 + j) * NSSY;

      if(m < 0 || m >= NSYM) continue;

      vmax = -1;
      nmax = -1;

      for(n = 0; n < 8; ++n)
      {
        v = map[m * NFFT + cand->i + n * NFOS];
        if(v > vmax)
        {
          vmax = v;
          nmax = n;
        }
      }
      if(nmax == costas[j]) ++result;
    }
  }
  return result;
}

/*
 * Symbol stage: log-likelihood ratios from the tone amplitudes.
 */

static inline real_t max(real_t a, real_t b, real_t c, real_t d)
{
  real_t x, y;
  x = a > b ? a : b;
  y = c > d ? c : d;
  return x > y ? x : y;
}

static inline real_t tone_amp(int m)
{
  return (m >= 0 && m < NSYM * NFFT) ? map[m] : 0;
}

static inline real_t floored_log10(real_t v)
{
  return v > 0 ? log10f(v) : -10;
}

static void process(sync_t *cand)
{
  int i, j, k, l;
  real_t a[8], d, sum, avg, sig;

  for(i = 0; i < 2; ++i)
  {
    for(j = 0; j < 29; ++j)
    {
      k = cand->j * NFFT + cand->i + (i * 36 + j + 7) * NSSY * NFFT;

      for(l = 0; l < 8; ++l)
      {
        a[l] = tone_amp(k + graymap[l] * NFOS);
      }

      k = i * 87 + j * 3;

      llr[k + 0] = floored_log10(max(a[4], a[5], a[6], a[7])) - floored_log10(max(a[0], a[1], a[2], a[3]));
      llr[k + 1] = floored_log10(max(a[2], a[3], a[6], a[7])) - floored_log10(max(a[0], a[1], a[4], a[5]));
      llr[k + 2] = floored_log10(max(a[1], a[3], a[5], a[7])) - floored_log10(max(a[0], a[2], a[4], a[6]));
    }
  }

  sum = 0;
  for(i = 0; i < N; ++i)
  {
    sum += llr[i];
  }
  avg = sum / N;

  sum = 0;
  for(i = 0; i < N; ++i)
  {
    d = llr[i] - avg;
    sum += d * d;
  }
  sig = sqrtf(sum / (N - 1)) / 4;

  for(i = 0; i < N; ++i)
  {
    llr[i] /= sig;
  }
}

/*
 * CRC stage: verify the 14-bit CRC of the 91 information bits.
 */

static int check()
{
  static const uint8_t poly[15] = {1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1 };
  uint8_t data[96] = {0};
  int i, j;

  memcpy(data, message, 77);

  for(i = 0; i < 82; ++i)
  {
    if(data[i])
    {
      for(j = 0; j < 15; ++j)
      {
        data[i + j] ^= poly[j];
      }
    }
  }

  for(i = 0; i < 14; ++i)
  {
    if(data[82 + i] != message[77 + i]) return 0;
  }

  return 1;
}

/*
 * LDPC stage: soft-decision decoding of the 174/91 code.
 */

static inline real_t tanh_approx(real_t x)
{
  real_t x2 = x * x;

  return x > 3.64 ? 1 : x < -3.64 ? -1 : (945 + (105 + x2) * x2) * x / (945 + (420 + 15 * x2) * x2);
}

static inline real_t atanh_approx(real_t x)
{
  real_t x2 = x * x;

  return (945 - (735 - 64 * x2) * x2) * x / (945 - (1050 - 225 * x2) * x2);
}

static int decode(int iterations)
{
  int i, j, k, l, iter, ibj, current, previous, counter;
  real_t x, tnm, tov[N][3], toc[M][7], pre[8], suf[8];

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
      message[i] = llr[i] + tov[i][0] + tov[i][1] + tov[i][2] > 0;
    }

    current = 0;
    for(i = 0; i < M; ++i)
    {
      l = 0;
      for(j = 0; j < nrw[i]; ++j) l += message[nm[i][j]];
      if(l & 1) ++current;
    }

    if(current == 0 && check()) return 1;

    if(iter > 0)
    {
      counter = current < previous ? 0 : counter + 1;
      if(counter > 4 && iter > 9 && current > 15) return 0;
    }

    previous = current;

    for(ibj = 0; ibj < N; ++ibj)
    {
      tnm = llr[ibj] + tov[ibj][0] + tov[ibj][1] + tov[ibj][2];
      for(k = 0; k < 3; ++k)
      {
        x = (tov[ibj][k] - tnm) * 0.5;
        toc[mn[ibj][k]][cpos[ibj][k]] = tanh_approx(x);
      }
    }

    for(i = 0; i < M; ++i)
    {
      pre[0] = 1;
      for(j = 0; j < nrw[i]; ++j) pre[j + 1] = pre[j] * toc[i][j];
      suf[nrw[i]] = 1;
      for(j = nrw[i] - 1; j >= 0; --j) suf[j] = suf[j + 1] * toc[i][j];
      for(j = 0; j < nrw[i]; ++j)
      {
        x = -pre[j] * suf[j + 1];
        tov[nm[i][j]][p2[i][j]] = 2 * atanh_approx(x);
      }
    }
  }

  return 0;
}

/*
 * Message stage: unpack the 91 information bits.
 */

static void trim(char *s)
{
  char *p = s, *e = s + strlen(s);

  while(e > s && isspace(e[-1])) --e;
  *e = 0;
  while(p < e && isspace(*p)) ++p;
  memmove(s, p, e - p + 1);
}

/* Value of n consecutive message bits ending at message[top] (LSB first). */
static uint64_t bits(int top, int n)
{
  uint64_t v = 0;
  int i;

  for(i = 0; i < n; ++i) v |= (uint64_t)message[top - i] << i;
  return v;
}

static int unpack(char *call, char *grid)
{
  int i, n;
  uint64_t icall;
  uint16_t igrid;
  uint8_t i3;

  call[0] = 0;
  grid[0] = 0;

  i3 = bits(76, 3);

  if(i3 == 1 || i3 == 2)
  {
    icall = bits(56, 28);

    igrid = bits(73, 15);

    if(igrid <= MAXGRID && igrid != 32373)
    {
      n = igrid;
      grid[4] = 0;
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
    if(message[70]) return 0;

    icall = bits(69, 58);
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

static int snr(sync_t *cand)
{
  return floor(20.0 * log10f(1e-32 + cand->r) - 26 + 0.5);
}

int main(int argc, char **argv)
{
  FILE *fp;
  double dialfreq;
  int i, j, k, freq;
  sync_t *curr, *next, temp;
  char *date, *time, *suffix, call[12], grid[5];
  real_t dt, a[4] = {0.35875, 0.48829, 0.14128, 0.01168};

  if(argc != 2)
  {
    return EXIT_FAILURE;
  }

  if((fp = fopen(argv[1], "rb")) == NULL)
  {
    fprintf(stderr, "Cannot open input file %s.\n", argv[1]);
    return EXIT_FAILURE;
  }

  suffix = strstr(argv[1], ".c2");
  *suffix = 0;
  time = suffix - 4;
  *(suffix - 5) = 0;
  date = suffix - 11;

  signal = malloc(sizeof(complex_t) * NSAMP);
  map = malloc(sizeof(real_t) * NSYM * NFFT);
  list = malloc(sizeof(sync_t) * NFFT);

  for(i = 0; i < M; ++i)
  {
    nrw[i] = nm[i][6] < 0 ? 6 : 7;

    for(j = 0; j < nrw[i]; ++j)
    {
      for(k = 0; k < 3; ++k)
      {
        if(mn[nm[i][j]][k] == i)
        {
          p2[i][j] = k;
          break;
        }
      }
    }
  }

  for(i = 0; i < N; ++i)
  {
    for(k = 0; k < 3; ++k)
    {
      for(j = 0; j < nrw[mn[i][k]]; ++j)
      {
        if(nm[mn[i][k]][j] == i)
        {
          cpos[i][k] = j;
          break;
        }
      }
    }
  }

  buffer = pffft_aligned_malloc(sizeof(complex_t) * NFFT);
  setup = pffft_new_setup(NFFT, PFFFT_COMPLEX);

  for(i = 0; i < NSPS; ++i)
  {
    window[i] = a[0] -
      a[1] * cosf(2.0 * M_PI * i / (NSPS - 1)) +
      a[2] * cosf(4.0 * M_PI * i / (NSPS - 1)) -
      a[3] * cosf(6.0 * M_PI * i / (NSPS - 1));
  }

  fread(&dialfreq, 1, 8, fp);

  for(i = 0; i < 4; ++i)
  {
    fread(signal, 1, sizeof(complex_t) * NSAMP, fp);

    sync();

    for(j = 2; j < NFFT - 2; ++j)
    {
      curr = &list[j];
      next = &list[j + 1];

      if(curr->k == 0 || curr->s < SYNC_MIN) continue;

      if(next->k != 0 && next->s > curr->s)
      {
        temp = *curr;
        *curr = *next;
        *next = temp;
      }

      if(count_costas(curr) < 6) continue;

      process(curr);

      if(!decode(curr->s < SYNC_ITER_CAP ? LDPC_ITER_LOW : LDPC_ITER_FULL) || !unpack(call, grid)) continue;

      next->k = 0;

      dt = curr->j * NSTP / 4.0e3 - 0.5;
      freq = floor(dialfreq + curr->i * 4.0e3 / NFFT - 2.0e3 + 0.5);
      printf("%6s %4s%02d %5.1f %3d %5.2f %8d %11s %4s\n", date, time, i * 15, curr->s, snr(curr), dt, freq, call, grid);
    }
  }

  return EXIT_SUCCESS;
}
