#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// variáveis
static uint32_t seed;
static double  *A, *R, *Q;
static int      M, N;
static int      tamanho_sub_arranjo;
static int      tamanho_global, rank_global;
static int      linhas_processo;
static double  *vetor;

// Funções
static void gramschmidt_mpi(const int *inicio, const int *fim);
static void ler_args(int *argc, char ***argv);
static void preencher_matriz(void);
static void print_ajuda(void);
static void atualizar_matrizes(void);
static void escrever_matriz(void);

int
main(int argc, char **argv) {
  // Inicializações
  ler_args(&argc, &argv);
  preencher_matriz();

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &tamanho_global);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);

  linhas_processo = M / tamanho_global;
  int inicio      = rank_global * linhas_processo;
  int fim = (rank_global == tamanho_global - 1) ? M : inicio + linhas_processo;

  tamanho_sub_arranjo = linhas_processo * N;

  vetor = (double *) malloc(sizeof(double) * tamanho_sub_arranjo);

  gramschmidt_mpi(&inicio, &fim);

  MPI_Finalize();

  if (rank_global == 0) {
    escrever_matriz();
  }

  return 0;
}

static void
gramschmidt_mpi(const int *inicio, const int *fim) {
  double nrm, raiz_nrm;
  int    i, j, k, offsetA, offsetR, offsetQ;

  for (k = 0; k < N; k++) {
    nrm      = 0;
    raiz_nrm = 0;

    for (i = *inicio; i < *fim; i++) {
      offsetA  = i + M * k;
      nrm     += A[offsetA] * A[offsetA];
      if (rank_global > 0) {
        double nrm_aux;
        MPI_Recv(&nrm_aux, 1, MPI_DOUBLE, rank_global - 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        nrm += nrm_aux;
      }
      if (rank_global < tamanho_global - 1) {
        MPI_Send(&nrm, 1, MPI_DOUBLE, rank_global + 1, 0, MPI_COMM_WORLD);
      }
    }

    if (rank_global == tamanho_global - 1) {
      raiz_nrm = sqrt(nrm);
    }

    MPI_Bcast(&raiz_nrm, 1, MPI_DOUBLE, tamanho_global - 1, MPI_COMM_WORLD);
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
    if (tamanho_global != 1) {
      atualizar_matrizes();
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

static void
atualizar_matrizes(void) {
  int offset;
  if (rank_global < tamanho_global - 1) {
    offset = rank_global * linhas_processo * N;
    memcpy(vetor, A + offset, tamanho_sub_arranjo);
    MPI_Send(vetor, tamanho_sub_arranjo, MPI_DOUBLE, tamanho_global - 1, 0,
             MPI_COMM_WORLD);
  } else if (rank_global == tamanho_global - 1) {
    for (int i = 0; i < rank_global; i++) {
      offset = i * linhas_processo * N;
      MPI_Recv(vetor, tamanho_sub_arranjo, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      memcpy(A + offset, vetor, tamanho_sub_arranjo);
    }
  }

  MPI_Bcast(A, M * N, MPI_DOUBLE, tamanho_global - 1, MPI_COMM_WORLD);
}

static void
escrever_matriz(void) {
  FILE *arq = fopen("saida.txt", "w");
  if (arq == NULL) {
    perror("Erro ao abrir arquivo de saída");
    exit(EXIT_FAILURE);
  }

  fprintf(arq, "%d %d\n", M, N);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N - 1; j++) {
      fprintf(arq, "%.2f ", A[i + M * j]);
    }
    fprintf(arq, "%.2f\n", A[i + M * (N - 1)]);
  }

  fclose(arq);
}
