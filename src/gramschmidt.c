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
  int seed;
  char tamanho;
} Programa;

/* Array initialization. */
static void init_array(int m, int n, DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                       DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                       DATA_TYPE POLYBENCH_2D(Q, M, N, m, n));

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int m, int n, DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                        DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                        DATA_TYPE POLYBENCH_2D(Q, M, N, m, n));

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static void kernel_gramschmidt(int m, int n,
                               DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                               DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                               DATA_TYPE POLYBENCH_2D(Q, M, N, m, n));

static void definir_programa(Programa *programa, int argc, char *argv[]);
static void print_ajuda(void);
static void *gramschmidt_par(void *arg);
static void chamar_gs_seq(int argc, char *argv[]);
static void chamar_gs_par(int num_threads, int argc, char *argv[]);

static pthread_spinlock_t sl;
static pthread_barrier_t barreira;
static pthread_barrier_t barreira2;

int main(int argc, char **argv) {
  Programa programa = {0};
  definir_programa(&programa, argc, argv);

  if (programa.num_threads == 1) {
    chamar_gs_seq(argc, argv);
  } else {
    chamar_gs_par(programa.num_threads, argc, argv);
  }

  return 0;
}

static double nrm_real;

static void *gramschmidt_par(void *arg) {
  const Args *args = (Args *)arg;
  int k, i, j;

  double nrm;

  double R[N][N], Q[M][N], A[M][N];

  for (k = 0; k < N; k++) { // Lê todos os vetores
    nrm = SCALAR_VAL(0.0);
    nrm_real = 0.0;

    for (i = args->inicio; i < args->fim; i++) // PARALELIZA AQUI
      nrm += A[i][k] *
             A[i][k]; // Normal é a soma dos quadrados de cada elemento do vetor

    pthread_barrier_wait(&barreira);

    pthread_spin_lock(&sl);
    nrm_real += nrm;              // BARREIRA E LOCK
    R[k][k] = SQRT_FUN(nrm_real); // R[k][k] guarda a raiz da normal
    pthread_spin_unlock(&sl);

    for (i = args->inicio; i < args->fim; i++) // PARALELIZA AQUI
      Q[i][k] = A[i][k] / R[k][k];

    for (j = k + 1; j < N; j++) {
      R[k][j] = SCALAR_VAL(0.0);

      for (i = args->inicio; i < args->fim; i++) // PARALELIZA AQUI
        R[k][j] += Q[i][k] * A[i][j];

      for (i = args->inicio; i < args->fim; i++) // PARELILZA AQUI
        A[i][j] = A[i][j] - Q[i][k] * R[k][j];
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
        programa->tamanho = 's';
      } else if (strncmp(optarg, "medium", 6) == 0) {
        programa->tamanho = 'm';
      } else if (strncmp(optarg, "large", 5) == 0) {
        programa->tamanho = 'l';
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

static void chamar_gs_seq(int argc, char *argv[]) {
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, M, N, m, n);
  POLYBENCH_2D_ARRAY_DECL(R, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Q, DATA_TYPE, M, N, m, n);

  /* Initialize array(s). */
  init_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R),
                     POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A),
                                    POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);
}

static void chamar_gs_par(int num_threads, int argc, char *argv[]) {
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

static void init_array(int m, int n, DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                       DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                       DATA_TYPE POLYBENCH_2D(Q, M, N, m, n)) {
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE)((i * j) % m) / m) * 100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}

static void print_array(int m, int n, DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                        DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                        DATA_TYPE POLYBENCH_2D(Q, M, N, m, n)) {
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
    }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
    }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}

static void kernel_gramschmidt(int m, int n,
                               DATA_TYPE POLYBENCH_2D(A, M, N, m, n),
                               DATA_TYPE POLYBENCH_2D(R, N, N, n, n),
                               DATA_TYPE POLYBENCH_2D(Q, M, N, m, n)) {
  int i, j, k;

  DATA_TYPE nrm;

#pragma scop
  for (k = 0; k < _PB_N; k++) {
    nrm = SCALAR_VAL(0.0);
    for (i = 0; i < _PB_M; i++)
      nrm += A[i][k] * A[i][k];
    R[k][k] = SQRT_FUN(nrm);
    for (i = 0; i < _PB_M; i++)
      Q[i][k] = A[i][k] / R[k][k];
    for (j = k + 1; j < _PB_N; j++) {
      R[k][j] = SCALAR_VAL(0.0);
      for (i = 0; i < _PB_M; i++)
        R[k][j] += Q[i][k] * A[i][j];
      for (i = 0; i < _PB_M; i++)
        A[i][j] = A[i][j] - Q[i][k] * R[k][j];
    }
  }
#pragma endscop
}
