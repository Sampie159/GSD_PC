/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

/* Include polybench common header. */
#include "polybench.h"

/* Include benchmark-specific header. */
#include "gramschmidt.h"

typedef struct {
  int inicio;
  int fim;
} Args;

typedef struct {
  int num_threads;
  uint32_t seed;
  char tamanho;
} Programa;

static int M, N;

/* Array initialization. */
static void init_array(uint32_t seed, int m, int n);

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int m, int n, double *A, double *R, double *Q);

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static void kernel_gramschmidt(int m, int n, double *A, double *R, double *Q);

static void definir_programa(Programa *programa, int argc, char *argv[]);
static void print_ajuda(void);
static void *gramschmidt_par(void *arg);
static void chamar_gs_par(int num_threads);

static pthread_spinlock_t sl;
static pthread_barrier_t barreira;
static pthread_barrier_t barreira2;

static double *A, *R, *Q;

int main(int argc, char **argv) {
  Programa programa = {0};
  definir_programa(&programa, argc, argv);

  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Initialize array(s). */
  init_array(programa.seed, m, n);

  /* Start timer. */
  /* polybench_start_instruments; */

  if (programa.num_threads == 1) {
    /* Run kernel. */
    kernel_gramschmidt(m, n, A, R, Q);
  } else {
    chamar_gs_par(programa.num_threads);
  }

  /* Stop and print timer. */
  /* polybench_stop_instruments; */
  /* polybench_print_instruments; */

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  /* polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), */
  /*                                   POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q))); */

  /* Be clean. */
  /* POLYBENCH_FREE_ARRAY(A); */
  /* POLYBENCH_FREE_ARRAY(R); */
  /* POLYBENCH_FREE_ARRAY(Q); */

  return 0;
}

static double nrm_real;

static void *gramschmidt_par(void *arg) {
  const Args *args = (Args *)arg;
  int k, i, j, offsetA, offsetR, offsetQ;

  double nrm;

  for (k = 0; k < N; k++) { // Lê todos os vetores
    nrm = SCALAR_VAL(0.0);
    nrm_real = 0.0;

    for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
      offsetA = i + M * k;
      nrm += A[offsetA] *
             A[offsetA]; // Normal é a soma dos quadrados de cada elemento do vetor
    } 

    pthread_spin_lock(&sl);
    nrm_real += nrm;              // BARREIRA E LOCK
    pthread_spin_unlock(&sl);

    offsetR = k + N * k;
    R[offsetR] = SQRT_FUN(nrm_real); // R[k][k] guarda a raiz da normal
    pthread_barrier_wait(&barreira);
    for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
      offsetA = i + M * k;
      offsetR = k + N * k;
      Q[offsetA] = A[offsetA] / R[offsetR];
    }

    for (j = k + 1; j < N; j++) {
      offsetR = k + N * j;
      R[offsetR] = SCALAR_VAL(0.0);

      for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
        offsetA = i + M * j;
        offsetQ = i + M * k;
        offsetR = k + N * j;
        R[offsetR] += Q[offsetQ] * A[offsetA];
      }

      for (i = args->inicio; i < args->fim; i++) { // PARELILZA AQUI
        offsetA = i + M * j;
        offsetQ = i + M * k;
        offsetR = k + N * j;
        A[offsetA] = A[offsetA] - Q[offsetQ] * R[offsetR];
      }
    }
    pthread_barrier_wait(&barreira2);
  }

  pthread_exit(NULL);
}

static void definir_programa(Programa *programa, int argc, char *argv[]) {
  int opt, aux;
  while ((opt = getopt(argc, argv, "d:t:S:h")) != -1) {
    switch (opt) {
    case 'd':
      if (strncmp(optarg, "small", 5) == 0) {
        M = 4000;
        N = 4800;
      } else if (strncmp(optarg, "medium", 6) == 0) {
        M = 5000;
        N = 6000;
      } else if (strncmp(optarg, "large", 5) == 0) {
        M = 6000;
        N = 7200;
      } else {
        printf("Tamanho não suportado, use: small, medium ou large\n");
        exit(EXIT_FAILURE);
      }
      break;

    case 't':
      aux = atoi(optarg);
      if (aux < 1) {
        printf("Número de threads tem que ser maior ou igual a 1\n");
        exit(EXIT_FAILURE);
      } else {
        programa->num_threads = aux;
      }
      break;
      
    case 'S':
      programa->seed = atoi(optarg);
      break;

    case 'h':
      print_ajuda();
      exit(EXIT_SUCCESS);
    }
  }
}

static void print_ajuda(void) {
  printf(
      "Uso: ./gsd -d <TAMANHO> -t <NUM_THREADS> -S <SEED>\n"
      "Opções:\n"
      "-d: <TAMANHO> small | medium | large - Define qual tamanho de matriz "
      "irá utilizar.\n"
      "-t: <NUM_THREADS>                    - Define a quantidade de threads "
      "que irá utilizar.\n"
      "-S: <SEED>                           - Define a seed que irá utilizar "
      "para gerar os valores das matrizes.\n"
      "-h:                                  - Escreve essa mensagem de "
      "ajuda.\n");
}

static void chamar_gs_par(int num_threads) {
  pthread_barrier_init(&barreira, NULL, num_threads);
  pthread_barrier_init(&barreira2, NULL, num_threads);
  pthread_spin_init(&sl, num_threads);

  pthread_t threads[num_threads];
  int qt_linha_thread = M / num_threads;

  for (int i = 0; i < num_threads; i++) {
    Args *args = (Args *)malloc(sizeof(Args));
    args->inicio = qt_linha_thread * i;
    args->fim = (i == num_threads - 1) ? M : args->inicio + qt_linha_thread;
    pthread_create(&threads[i], NULL, gramschmidt_par, (void *)args);
  }

  for (int i = 0; i < num_threads; i++) {
    pthread_join(threads[i], NULL);
  }

  pthread_barrier_destroy(&barreira);
  pthread_barrier_destroy(&barreira2);
  pthread_spin_destroy(&sl);
}

static void init_array(uint32_t seed, int m, int n) {
  A = (double *)malloc(M * N * sizeof(double));
  R = (double *)malloc(N * N * sizeof(double));
  Q = (double *)malloc(M * N * sizeof(double));
  int i, j, offset;
  srand(seed);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      offset = i + M * j;
      A[offset] = ((float)rand()/(float)RAND_MAX) * 1000;
      Q[offset] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      offset = i + M * j;
      R[offset] = 0.0;
}

/* static void print_array(int m, int n, double *A, double *R, double *Q) { */
  /* int i, j; */

  /* /1* POLYBENCH_DUMP_START; *1/ */
  /* /1* POLYBENCH_DUMP_BEGIN("R"); *1/ */
  /* for (i = 0; i < n; i++) */
  /*   for (j = 0; j < n; j++) { */
  /*     if ((i * n + j) % 20 == 0) */
  /*       fprintf(POLYBENCH_DUMP_TARGET, "\n"); */
  /*     fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]); */
  /*   } */
  /* /1* POLYBENCH_DUMP_END("R"); *1/ */

  /* /1* POLYBENCH_DUMP_BEGIN("Q"); *1/ */
  /* for (i = 0; i < m; i++) */
  /*   for (j = 0; j < n; j++) { */
  /*     if ((i * n + j) % 20 == 0) */
  /*       fprintf(POLYBENCH_DUMP_TARGET, "\n"); */
  /*     fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]); */
  /*   } */
  /* /1* POLYBENCH_DUMP_END("Q"); *1/ */
  /* /1* POLYBENCH_DUMP_FINISH; *1/ */
/* } */

static void kernel_gramschmidt(int m, int n, double *A, double *R, double *Q) {
  int i, j, k, offsetA, offsetR, offsetQ;

  DATA_TYPE nrm;

#pragma scop
  for (k = 0; k < _PB_N; k++) {
    nrm = SCALAR_VAL(0.0);
    for (i = 0; i < _PB_M; i++) {
      offsetA = i + M * k;
      nrm += A[offsetA] * A[offsetA];
    }
    R[k + M * k] = SQRT_FUN(nrm);
    for (i = 0; i < _PB_M; i++) {
      offsetA = i + M * k;
      offsetR = k + M * k;
      Q[offsetA] = A[offsetA] / R[offsetR];
    }
    for (j = k + 1; j < _PB_N; j++) {
      offsetR = k + M * j;
      R[offsetR] = SCALAR_VAL(0.0);
      for (i = 0; i < _PB_M; i++) {
        offsetA = i + M * j;
        offsetQ = i + M * k;
        offsetR = k + M * j;
        R[offsetR] += Q[offsetQ] * A[offsetA];
      }
      for (i = 0; i < _PB_M; i++) {
        offsetA = i + M * j;
        offsetQ = i + M * k;
        offsetR = k + M * j;
        A[offsetA] = A[offsetA] - Q[offsetQ] * R[offsetR];
      }
    }
  }
#pragma endscop
}
