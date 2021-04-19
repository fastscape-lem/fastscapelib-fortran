.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .e .r .y .yr .ye .l .p .sh .csh .h

ifneq ($(findstring geodyn,$(HOSTNAME)),)
  F90 = mpif90
  F90_OPTIONS = -O -Mfreeform -Mdclchk -Mstandard -Mpreprocess
else
  F90 = gfortran
  CC = gcc-10
  F90_OPTIONS = -ffree-form -ffree-line-length-none -fimplicit-none -fbounds-check -Wall -Wno-unused-dummy-argument -cpp -fall-intrinsics -fallow-argument-mismatch  # needed to remove -std=f2008 -ffree-line-length-none is required for long macros
endif

ON_WINDOWS = f

FASTSCAPE_OBJ = FastScapeErrorCodes.o FastScape_ctx.o FastScape_api.o StreamPowerLaw.o Diffusion.o LocalMinima.o Advect.o Marine.o VTK.o FlowRouting.o Strati.o Uplift.o TerrainDerivatives.o

.f90.o:
	$(F90) -c $(F90_OPTIONS) $*.f90

all: $(FASTSCAPE_OBJ) LIB
	make LIB
	make RUNF90
	make RUNC

LIB: $(FASTSCAPE_OBJ)
	ar rc libFastScape.a $(FASTSCAPE_OBJ)

RUNF90: LIB
	$(F90) -O FastScape_run.f90 -o FastScape_run_f90 -L. -lFastScape

RUNC: LIB
	$(CC) -O FastScape_run.c -o FastScape_run_c -L. -lFastScape -lgfortran

clean:
	rm *.o *.mod *.a
