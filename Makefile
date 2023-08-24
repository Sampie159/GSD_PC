CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -g -std=gnu17 -O3

SRCS = $(wildcard $(SRCDIR)/*.c)
SRCDIR = src
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))
OBJDIR = obj

all: gsd gsd_mpi

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

gsd: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

gsd_mpi: src_mpi/main.c
	mpicc $(CFLAGS) $^ -o $@

clean:
	rm -rf $(OBJDIR) gsd
