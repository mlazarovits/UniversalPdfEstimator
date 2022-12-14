#include ROOT cflags and libraries
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)

#set c(xx)flags and libraries
CXXFLAGS    = $(ROOTCFLAGS)

GLIBS       = $(ROOTGLIBS)
#specify compiler
CXX         = g++

#specify local paths
INCLUDEDIR  = ./include/
SRCDIR      = ./src/
#make sure compiler knows where local include files are
CXX	   += -I$(INCLUDEDIR) -I.
OUTOBJ	    = ./obj/

#specify local source, include, and object files
CC_FILES    = $(wildcard src/*.cc)
HH_FILES    = $(wildcard include/*.hh)
OBJ_FILES   = $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))


#specify what to make
all: toyBezier.x toyWeightedBezier.x toyBezierTransform.x
test: TEST_newToyBezier.x

#executables
toyBezier.x: $(SRCDIR)toyBezier.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o toyBezier.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch toyBezier.x

toyWeightedBezier.x: $(SRCDIR)toyWeightedBezier.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o toyWeightedBezier.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch toyWeightedBezier.x

toyBezierTransform.x: $(SRCDIR)toyBezierTransform.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o toyBezierTransform.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch toyBezierTransform.x

TEST_newToyBezier.x: $(SRCDIR)TEST_newToyBezier.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o TEST_newToyBezier.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch TEST_newToyBezier.x

#where to put object files
$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

#clean options
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -f AutoDict*
	rm -f *.d



