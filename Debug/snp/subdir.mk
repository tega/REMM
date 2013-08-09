################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../snp/snp1.o \
../snp/snp2.o \
../snp/snp3.o \
../snp/snp4.o \
../snp/snp5.o 

F90_SRCS += \
../snp/remm.f90 \
../snp/snp1.f90 \
../snp/snp2.f90 \
../snp/snp3.f90 \
../snp/snp4.f90 \
../snp/snp5.f90 

OBJS += \
./snp/remm.o \
./snp/snp1.o \
./snp/snp2.o \
./snp/snp3.o \
./snp/snp4.o \
./snp/snp5.o 


# Each subdirectory must supply rules for building sources it contributes
snp/%.o: ../snp/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

snp/remm.o: ../snp/remm.f90 snp/datatypes.o snp/startup.o snp/steps.o

snp/snp1.o: ../snp/snp1.f90

snp/snp2.o: ../snp/snp2.f90 snp/datatypes.o

snp/snp3.o: ../snp/snp3.f90 snp/datatypes.o snp/simulation.o

snp/snp4.o: ../snp/snp4.f90 snp/datatypes.o snp/simulation.o snp/likelihood_utilities.o

snp/snp5.o: ../snp/snp5.f90 snp/datatypes.o snp/startup.o snp/steps.o


