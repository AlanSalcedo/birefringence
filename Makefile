# make file to compile and link a test program that uses the verbmenu library


CCFILE       = test.cc birefringence.cc

# Define filename suffixes
ObjSuf        = .o
SrcSuf        = .cc
IncSuf        = .hh
ExeSuf        = .exe
DllSuf        = .so
OutPutOpt     = -o

# Define the compile and link commands
CXX           = g++
CXXFLAGS      = -c -O3 -pg -o $(OBJ) -Wall -fPIC -I$(ROOTSYS)/include -I${GSLDIR}/include
#-I$(HOOVERINC) 
LD            = g++
LDFLAGS       = -O3 -pg -L{GSLDIR}/lib

# Define root libraries
#ROOTLIBS      = -L$(ROOTSYS)/lib 
ROOTLIBS      = `${ROOTSYS}/bin/root-config --libs` 

# Define all my libraries	
LIBS          = $(ROOTLIBS) -lMinuit -lgsl -lgslcblas

# Define shortcuts for compiling and linking
COMPILE = $(CXX) $(CXXFLAGS) $(SCRATCH)
LINK = 	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) $(OutPutOpt) $(EXE)
#------------------------------------------------------------------------------

# Define the file names
OBJ       = test.obj
SRC       = $(CCFILE)
EXE       = test.exe
INC       = test.inc 
SCRATCH	  = tmp.cc
#HOOVERINC    = /home/hoover/icemc/

all:            $(EXE)

# Update the executable if the object file has changed
$(EXE): 	$(OBJ)
		$(LINK)


# Update the object file if the source, or include file changed
$(OBJ):     $(SRC)
		cat $(CCFILE) > $(SCRATCH)
		$(COMPILE)

clean:
	@rm -f $(OBJ) $(SCRATCH) core

.SUFFIXES: .$(SrcSuf)

.$(SrcSuf).$(ObjSuf):	
	$(CXX) $(CXXFLAGS) -c $<	







