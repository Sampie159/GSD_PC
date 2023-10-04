CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -g -std=gnu17 -O3 -pipe -march=native

SRCS = $(wildcard $(SRCDIR)/*.c)
SRCDIR = src
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))
OBJDIR = obj

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@ -fopenmp

$(OBJDIR):
	mkdir -p $(OBJDIR) 

gsd: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm -fopenmp

gsd_mpi: src_mpi/main.c
	mpicc $(CFLAGS) $^ -o $@ -lm

gsd_ompi: src_mpi/main.c
	mpicc $(CFLAGS) $^ -o $@ -lm -D_OPENMPI

clean:
	rm -rf $(OBJDIR) gsd gsd_mpi
