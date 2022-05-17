INCLUDES = -I "./include"
CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -mfpmath=387
OBJ = tree_based.o

tree_based:
	$(CC) $(CFLAGS) math_funcs.c tree_based_heston.c $(LIBS) -o ./bin/$@

tree_based_ibm:
	$(CC) $(CFLAGS) math_funcs.c tree_based_ibm.c $(LIBS) -o ./bin/$@

tree_based_sp:
	$(CC) $(CFLAGS) math_funcs.c tree_based_sp.c $(LIBS) -o ./bin/$@

quad_tree:
	$(CC) $(CFLAGS) $(LIBS) math_funcs.c quadrinomial_tree.c -o ./bin/$@

quad_tree_amer:
	$(CC) $(CFLAGS) $(LIBS) math_funcs.c quad_tree_amer.c -o ./bin/$@

quad_tree_inter:
	$(CC) $(CFLAGS) $(LIBS) math_funcs.c quad_tree_interlace.c -o ./bin/$@

quad_tree_sp:
	$(CC) $(CFLAGS) $(LIBS) math_funcs.c quad_tree_sp.c -o ./bin/$@