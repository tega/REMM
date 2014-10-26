#Define the fortran 77/90 compiler
FF = gfortran 

#Define the NPSOL library path
NPSOL_LIB = ../../npsol/lib

#Define the objects files
OBJS := snp1.o snp2.o snp3.o snp4.o SNP.o SNP_ROBUST.o remm.o

remm: $(OBJS) 
	$(FF) -o remm $(OBJS) $(NPSOL_LIB)/libnpsol.a
	mv remm ../run/remm

snp1.o : snp1.f90
	$(FF) -Wall -c snp1.f90
snp2.o : snp2.f90 snp1.o  
	$(FF) -Wall -c snp2.f90
snp3.o : snp3.f90 snp2.o snp1.o 
	$(FF) -Wall -c snp3.f90
snp4.o : snp4.f90 snp3.o snp2.o snp1.o
	$(FF) -Wall -c snp4.f90
SNP.o : SNP.f90 snp4.o snp3.o snp2.o snp1.o
	$(FF) -Wall -c SNP.f90
SNP_ROBUST.o : SNP_ROBUST.f90 snp4.o snp3.o snp2.o snp1.o
	$(FF) -Wall -c SNP_ROBUST.f90
remm.o : remm.f90 SNP_ROBUST.o SNP.o snp4.o snp3.o snp2.o snp1.o
	$(FF) -Wall -c remm.f90


clean:
	rm $(OBJS) *.mod
	rm ../run/remm

