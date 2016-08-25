#
# Makefile for MC4LJsystem
#

TARGET = MD4LJsystem

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = MD_LJ.f90 force_LJ.f90 neibhoringlist.f90 ini_conf.f90 gen_velo.f90 mtmod.f90 \
	 pbs.f90 pdb.f90 read_inputs.f90 write_outputs.f90 reallocate.f90

OBJS = MD_LJ.o force_LJ.o neibhoringlist.o ini_conf.o gen_velo.o mtmod.o \
	pbs.o pdb.o read_inputs.o write_outputs.o reallocate.o

MOD_FILES = force_LJ.mod neibhoringlist.mod ini_conf.mod gen_velo.mod mtmod.mod \
	 pbs.mod pdb.mod read_inputs.mod write_outputs.mod reallocate.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

MD_LJ.o: ini_conf.o mtmod.o gen_velo.o pbs.o force_LJ.o reallocate.o neibhoringlist.o pdb.o read_inputs.o write_outputs.o
force_LJ.o:
neiboringlist.o: reallocate.o
ini_conf.o:
gen_velo.o: mtmod.o
mtmod.o:
pbs.o:
pdb.o:
read_input.o: pdb.o
write_output.o: pdb.o
reallocate.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
