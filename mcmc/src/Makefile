CLIB = -lm

CC = gcc

OBJS1 = chain_covar.o readpar.o hod.o getobs.o random.o mcmc_util.o matrix.o file_util.o table_util.o nrutil.o mathutil.o 
chain_covar: $(OBJS1)
	$(CC) -o chain_covar $(OBJS1) $(CLIB)

OBJS2 = chain_covar.o readpar.o hod.o getobs_cscorr.o random.o mcmc_util.o matrix.o file_util.o table_util.o nrutil.o mathutil.o
chain_covar_cscorr: $(OBJS2)
	$(CC) -o chain_covar_cscorr $(OBJS2) $(CLIB)


clean:
	rm -f *.o

