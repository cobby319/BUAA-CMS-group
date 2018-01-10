EXEC     = makehistos\
	#Add other executables here if needed...
CXX      = g++
CXXFLAGS = -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include
SRC      = $(wildcard Sources/*.cc)
OBJ      = $(SRC:.cc=.o)

all: $(EXEC)

makehistos: makehistos.cc $(OBJ)
	$(CXX) $< $(CXXFLAGS) $(OBJ) -o $@

%.o: %.cc %.h
	$(CXX) -c $< $(CXXFLAGS) -o $@

.PHONY: clean mrproper

clean:
	rm -rf *.o *.d Sources/*.o

mrproper: clean
	rm -rf $(EXEC)
