

include make.inc



# CHANGEME: This should be the name of your executable
EXE = iquad
LIBNAME = libiquad.a
MYLIBSDIR = lib
MYINCSDIR = include


# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS_LIB = tools.o biqmac.o boxqp.o qpgen.o quadprog.o myqp.o sdpVectors.o newbb.o refProb2.o runProblem.o

OBJS = $(OBJS_LIB) ampl.o main.o


OBJS_SERVER = 
EXE_SERVER = 

OBJS_CLIENT =
EXE_CLIENT = 


# Additional libraries
ADDLIBS = 

# Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =





# additional C++ Compiler options for linking

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)

INCL = $(CSDPINC) $(SPMINC) $(MINLPPROBINC) $(BBLINC) $(OPTSINC) $(ASLINC)   $(ADDINCFLAGS)



# Linker flags
#LIBS = $(CSDPLIB) $(MURIQUILIB) $(BBLLIB) $(OPTSUBLIB) $(SPMLIB) $(MINLPPROBSUBLIB)  $(LAPACKLIB)


ifeq ("$(CXX)" , "cl")
	
	LIBS = $(BBLLIBNAME)  $(OPTSLIBNAME) $(SPMLIBNAME)  $(MINLPPROBLIBNAME)
else
# g++ and icpc
	LIBS = $(BBLLIB) $(OPTSSUBLIB) $(SPMLIB)  $(MINLPPROBLIB)
endif






# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: lib cleanexe $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS) sublibs
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@  $(OBJS)  $(ADDLIBS) $(LIBS)

clean: cleanexe
	- $(RM_CMD) $(OBJS) $(LIBNAME)

cleanexe:
	- $(RM_CMD) $(EXE) $(EXE).exe  $(EXE_CLIENT) $(EXE_CLIENT).exe

cleanall: clean  cleansublibs
	cd $(SPMDIR); make clean;
	cd $(MINLPPROBDIR); make clean;
	cd $(BBLDIR); make clean;
	cd $(OPTSDIR); make clean;

lib: $(OBJS_LIB)
	$(AR_CMD)$(LIBNAME)  $(OBJS_LIB)

sublibs:
	$(MAKE_CMD) -C $(MINLPPROBDIR)
	$(MAKE_CMD) -C $(BBLDIR)
	$(MAKE_CMD) -C $(OPTSDIR)

cleansublibs:
	$(MAKE_CMD) clean -C $(MINLPPROBDIR)
	
	$(MAKE_CMD) clean -C  $(BBLDIR)
	
	$(MAKE_CMD) clean -C  $(OPTSDIR) 

install: 
	
	- $(MKDIR_CMD) $(MYLIBSDIR)
	- $(MKDIR_CMD) $(MYINCSDIR)
	
	- $(COPY_CMD) $(MINLPPROBDIR)$(SYS_SEP)libminlpprob.a $(MYLIBSDIR)$(SYS_SEP)libminlpprob.a;
	- $(COPY_CMD) $(MINLPPROBDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(BBLDIR)$(SYS_SEP)libbbl.a $(MYLIBSDIR)$(SYS_SEP)libbbl.a;
	- $(COPY_CMD) $(BBLDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(OPTSDIR)$(SYS_SEP)liboptsolvers.a $(MYLIBSDIR)$(SYS_SEP)liboptsolvers.a;
	- $(COPY_CMD) $(OPTSDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(SPMDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(LIBNAME) $(MYLIBSDIR)$(SYS_SEP)$(LIBNAME)
	- $(COPY_CMD) *.h* $(MYINCSDIR)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ '$<'
