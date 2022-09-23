WFADIR=include/WFA2-lib
HCLUST=include/hclust-cpp


INC = -L$(WFADIR)/lib -I$(WFADIR) -L${HCLUST} -I${HCLUST}
CFLAGS=-lz -lm -fopenmp -lhts -lwfacpp
LDFLAGS=$(CFLAGS)
CXX=g++
CXX_FLAGS=-std=c++17 -Wall -O2 -g

ODIR = build
SDIR = src


OUT = organo

SRCFILES := $(wildcard src/*.cpp)

OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCFILES))

DEPS = $(OBJS:%.o=%.d)

.PHONY: all clean

all: $(ODIR)/$(OUT)

clean:
	$(RM) $(ODIR)/$(OUT) $(OBJS) $(DEPS)

packages:
	cd include/; git clone https://github.com/smarco/WFA2-lib.git; cd WFA2-lib; git reset --hard 06af50ce1b7a0d477a439a7ac901e4a18d76397a; make clean setup lib_wfa

$(ODIR)/$(OUT): $(OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(INC) $(CFLAGS)

-include $(DEPS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CXX_FLAGS) -MMD -MP -c -o $@ $< $(INC) $(LDFLAGS)