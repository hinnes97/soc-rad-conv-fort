FC=gfortran
LD=gfortran
SRCDIR=src
FCFLAGS= -g -fbacktrace -fcheck=all -Wall -O2
FLFLAGS= -L/usr/local/lib -llapack -lrefblas
OBJDIR=obj
OBJS=$(patsubst $(SRCDIR)/%.f90, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.f90))

$(info $$OBJS is [${OBJS}])

.PHONY: clean all depend

all: main.exe

main.exe: $(OBJS)
	$(LD) -o $@ $^ `nf-config --flibs` $(FLFLAGS)

$(OBJDIR)/%.o $(OBJDIR)/%.mod: $(SRCDIR)/%.f90
	$(FC) -c -J$(OBJDIR) -I$(OBJDIR) $(FCFLAGS) `nf-config --fflags` $< -o $@ 

include .depend

clean:
	rm -rf $(OBJDIR)/* main.exe .depend

depend .depend:
	makedepf90 -b $(OBJDIR) $(SRCDIR)/*.f90 > .depend
