#include <math.h>
#include <mpi/mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// TODO: Melhorar definição de Slice
typedef struct {
  int     tamanho;
  double *vetor;
} Slice;

// variáveis
static uint32_t seed;
static double  *A, *R, *Q;
static int      M, N;

// Funções
static void gramschmidt_mpi(const int *inicio, const int *fim,
                            const int *tamanho_global, const int *rank_global,
                            const int *linhas_processo);
static void ler_args(int *argc, char ***argv);
static void preencher_matriz(void);
static void print_ajuda(void);

int
main(int argc, char **argv) {
  // Inicializações
  ler_args(&argc, &argv);
  preencher_matriz();

  int tamanho_global, rank_global;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &tamanho_global);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);

  int linhas_processo = M / tamanho_global;
  int inicio          = rank_global * linhas_processo;
  int fim = (rank_global == tamanho_global - 1) ? M : inicio + linhas_processo;

  gramschmidt_mpi(&inicio, &fim, &tamanho_global, &rank_global,
                  &linhas_processo);

  MPI_Finalize();

  return 0;
}

static void
gramschmidt_mpi(const int *inicio, const int *fim, const int *tamanho_global,
                const int *rank_global, const int *linhas_processo) {
  double nrm, raiz_nrm;
  int    i, j, k, offsetA, offsetR, offsetQ;

  for (k = 0; k < N; k++) {
    nrm = 0;

    for (i = *inicio; i < *fim; i++) {
      offsetA  = i + M * k;
      nrm     += A[offsetA] * A[offsetA];
      if (rank_global > 0) {
        double nrm_aux;
        MPI_Recv(&nrm_aux, 1, MPI_DOUBLE, *rank_global - 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        nrm += nrm_aux;
      }
      if (rank_global < tamanho_global - 1) {
        MPI_Send(&nrm, 1, MPI_DOUBLE, *rank_global + 1, 0, MPI_COMM_WORLD);
      }
    }

    if (rank_global == tamanho_global - 1) {
      raiz_nrm = sqrt(nrm);
    }

    MPI_Bcast(&raiz_nrm, 1, MPI_DOUBLE, *tamanho_global - 1, MPI_COMM_WORLD);
    offsetR    = k + N * k;
    R[offsetR] = raiz_nrm;

    for (i = *inicio; i < *fim; i++) {
      offsetA    = i + M * k;
      Q[offsetA] = A[offsetA] / raiz_nrm;
    }

    for (j = k + 1; j < N; j++) {
      offsetR    = k + N * j;
      R[offsetR] = 0;

      for (i = *inicio; i < *fim; i++) {
        offsetA     = i + M * k;
        offsetQ     = i + M * j;
        R[offsetR] += Q[offsetA] * A[offsetQ];
      }

      for (i = *inicio; i < *fim; i++) {
        offsetA     = i + M * k;
        offsetQ     = i + M * j;
        A[offsetQ] -= R[offsetR] * Q[offsetA];
      }
    }

    if (rank_global > 0) {
      // TODO: Enviar slice de A para o próximo processo
    }
  }
}

static void
ler_args(int *argc, char ***argv) {
  int opt;
  while ((opt = getopt(*argc, *argv, "d:S:h")) != -1) {
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
        perror("Opção inválida de tamanho");
        exit(EXIT_FAILURE);
      }
      break;

    case 'S': seed = atoi(optarg); break;

    case 'h': print_ajuda(); break;
    }
  }
}

static void
preencher_matriz(void) {
  A = (double *) malloc(M * N * sizeof(double));
  R = (double *) malloc(N * N * sizeof(double));
  Q = (double *) malloc(M * N * sizeof(double));

  srand(seed);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      A[i + M * j] = (double) rand() / (double) RAND_MAX;
    }
  }
}

static void
print_ajuda(void) {
  printf(
      "Uso: mpirun -np <numero_processos> ./gsd_mpi [OPÇÕES]\n"
      "Opções:\n"
      "  -d <small|medium|large>  Tamanho da matriz (default: small)\n"
      "  -S <seed>                Semente para o gerador de números "
      "aleatórios\n"
      "  -h                       Exibe esta mensagem de ajuda\n");
}
