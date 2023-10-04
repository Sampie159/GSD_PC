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
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
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
  int      num_threads;
  uint32_t seed;
} Programa;

static int M, N;

/* Array initialization. */
static void init_array(uint32_t seed, int m, int n);

static void  definir_programa(Programa *programa, int argc, char *argv[]);
static void  print_ajuda(void);
static void *gramschmidt_par(void *arg);
static void  gram_schmidt(int num_threads);

static pthread_spinlock_t lock;
static pthread_barrier_t  barreira;
static pthread_barrier_t  barreira2;
static pthread_barrier_t  barreira3;
static pthread_barrier_t  barreira4;

static double *A, *R, *Q;

int
main(int argc, char **argv) {
  Programa programa = { 0 };
  definir_programa(&programa, argc, argv);

  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Initialize array(s). */
  init_array(programa.seed, m, n);

  double inicio = omp_get_wtime();
  gram_schmidt(programa.num_threads);
  double fim = omp_get_wtime();

  printf("Tempo: %f\n", fim - inicio);

  free(A);
  free(R);
  free(Q);

  return 0;
}

static double nrm;

static void *
gramschmidt_par(void *arg) {
  Args *args = (Args *) arg;
  int   k, i, j, offsetA, offsetR, offsetQ;

  printf("inicio: %d, fim: %d\n", args->inicio, args->fim);

  for (k = 0; k < N; k++) { // Lê todos os vetores
    nrm = SCALAR_VAL(0.0);

    pthread_barrier_wait(&barreira2);

    for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
      offsetA = i + M * k;
      pthread_spin_lock(&lock);
      nrm += A[offsetA] * A[offsetA]; // Normal é a soma dos quadrados de cada
                                      // elemento do vetor
      pthread_spin_unlock(&lock);
    }

    offsetR    = k + N * k;
    R[offsetR] = SQRT_FUN(nrm); // R[k][k] guarda a raiz da normal

    pthread_barrier_wait(&barreira);

    for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
      offsetA    = i + M * k;
      Q[offsetA] = A[offsetA] / R[offsetR];
    }

    for (j = k + 1; j < N; j++) {
      offsetR    = k + N * j;
      R[offsetR] = SCALAR_VAL(0.0);

      pthread_barrier_wait(&barreira4);

      for (i = args->inicio; i < args->fim; i++) { // PARALELIZA AQUI
        offsetA = i + M * j;
        offsetQ = i + M * k;
        pthread_spin_lock(&lock);
        R[offsetR] += Q[offsetQ] * A[offsetA];
        pthread_spin_unlock(&lock);
      }

      pthread_barrier_wait(&barreira3);

      for (i = args->inicio; i < args->fim; i++) { // PARELILZA AQUI
        offsetA     = i + M * j;
        offsetQ     = i + M * k;
        A[offsetA] -= Q[offsetQ] * R[offsetR];
      }
    }
  }

  free(args);

  return NULL;
}

static void
definir_programa(Programa *programa, int argc, char *argv[]) {
  int opt, aux;
  while ((opt = getopt(argc, argv, "d:t:S:h")) != -1) {
    switch (opt) {
    case 'd':
      if (strncmp(optarg, "small", 5) == 0) {
        M = 400;
        N = 480;
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

    case 'S': programa->seed = atoi(optarg); break;

    case 'h': print_ajuda(); exit(EXIT_SUCCESS);
    }
  }
}

static void
print_ajuda(void) {
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

static void
gram_schmidt(int num_threads) {
  pthread_barrier_init(&barreira, NULL, num_threads);
  pthread_barrier_init(&barreira2, NULL, num_threads);
  pthread_barrier_init(&barreira3, NULL, num_threads);
  pthread_barrier_init(&barreira4, NULL, num_threads);
  pthread_spin_init(&lock, 0);

  pthread_t threads[num_threads - 1];
  int       qt_linha_thread = M / num_threads;

  for (int i = 0; i < num_threads - 1; i++) {
    Args *args   = (Args *) malloc(sizeof(Args));
    args->inicio = qt_linha_thread * i;
    args->fim    = args->inicio + qt_linha_thread;
    pthread_create(&threads[i], NULL, gramschmidt_par, (void *) args);
  }

  Args *args   = (Args *) malloc(sizeof(Args));
  args->inicio = qt_linha_thread * (num_threads - 1);
  args->fim    = M;
  gramschmidt_par((void *) args);

  for (int i = 0; i < num_threads - 1; i++) {
    pthread_join(threads[i], NULL);
  }

  pthread_barrier_destroy(&barreira);
  pthread_barrier_destroy(&barreira2);
  pthread_barrier_destroy(&barreira3);
  pthread_barrier_destroy(&barreira4);
  pthread_spin_destroy(&lock);
}

static void
init_array(uint32_t seed, int m, int n) {
  A = (double *) malloc(M * N * sizeof(double));
  R = (double *) malloc(N * N * sizeof(double));
  Q = (double *) malloc(M * N * sizeof(double));
  int i, j, offset;
  srand(seed);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      offset    = i + M * j;
      A[offset] = ((float) rand() / (float) RAND_MAX) * 1000;
      Q[offset] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      offset    = i + M * j;
      R[offset] = 0.0;
    }
  }
}
