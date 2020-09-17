CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)
LIBS 		= $(ROOTLIBS)

ASSIMPLIBS       = $(ASSIMP)/lib/libassimp.so
ASSIMPINC        = $(ASSIMP)/include
LIBS 		+= -lm $(ASSIMPLIBS)
CXXFLAGS	+= -I$(ASSIMPINC)

CPPSRC = Global.cc Geometry.cc Material.cc Source.cc Triangle.cc Simulator.cc
CPPOBJ = $(CPPSRC:.cc=.o)

TARGET=OptSim
OBJS=OptSim.o

$(TARGET): $(OBJS) $(CPPOBJ) $(ASSIMPINC)
	@echo "Now making $@"
	@$(CPP) -o $@ $(OBJS) $(CPPOBJ) $(CXXFLAGS) $(LIBS)
	@echo "Compile done! \(^o^)/"

.cc.o:
	@echo "Compiling object file $<"
	@$(CPP) -c $(CXXFLAGS) $<

clean: 
	@echo "Now cleaning up"
	rm -f $(TARGET) *~ *.o *.o~ core
