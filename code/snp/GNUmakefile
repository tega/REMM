#Define the fortran 77/90 compiler
FF = gfortran -O3 

#Define the NPSOL directory
NPSOL_DIR = ../npsol/
LAPACK_DIR = ../lapack/

#Define the objects files
OBJS := snp1.o snp2.o snp3.o snp4.o snp5.o remm.o

remm : subsystem $(OBJS) 
	$(FF) -o remm $(OBJS) ../npsol/libnpsol.a ../lapack/libmylapack.a
	rm $(OBJS) *.mod
subsystem:
	$(MAKE) -f GNUmakefile_stat -C $(NPSOL_DIR)
	$(MAKE) -f GNUmakefile_stat -C $(LAPACK_DIR)
snp1.o : snp1.f90  
	$(FF) -Wall -c snp1.f90
snp2.o : snp2.f90 snp1.o  
	$(FF) -Wall -c snp2.f90
snp3.o : snp3.f90 snp2.o snp1.o 
	$(FF) -Wall -c snp3.f90
snp4.o : snp4.f90 snp3.o snp2.o snp1.o
	$(FF) -Wall -c snp4.f90
snp5.o : snp5.f90 snp4.o snp3.o snp2.o snp1.o
	$(FF) -Wall -c snp5.f90
remm.o : remm.f90 snp5.o snp4.o snp3.o snp2.o snp1.o
	$(FF) -Wall -c remm.f90
.PHONY : clean
clean: 
	rm remm
	$(MAKE) clean -f GNUmakefile_stat -C $(NPSOL_DIR)
	$(MAKE) clean -f GNUmakefile_stat -C $(LAPACK_DIR)