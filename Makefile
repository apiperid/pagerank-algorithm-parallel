SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc -O3
LM = -lm
FCILK = -fcilkplus
RM = rm -f

# ==========================================
# TARGETS

EXECUTABLES = PageRankSequential PageRankParallel test


all: $(EXECUTABLES)


PageRankSequential: PageRankSequential.c
	$(CC) $< -o $@ $(LM)

PageRankParallel: PageRankParallel.c
	$(CC) $< -o $@ $(LM) $(FCILK)

test: test.c
	$(CC) $< -o $@

clean:
	$(RM) *.o *~ $(EXECUTABLES)
