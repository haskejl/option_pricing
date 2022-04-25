INCLUDES = -I "./include"
CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -mfpmath=387
OBJ = tree_based.o

tree_based: tree_based_heston.c
	$(CC) $(CFLAGS) math_funcs.c tree_based_heston.c $(LIBS) -o ./bin/$@

quad_tree:
	$(CC) $(CFLAGS) $(LIBS) math_funcs.c quadrinomial_tree.c -o ./bin/$@

