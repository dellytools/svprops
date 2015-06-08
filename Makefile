DEBUG ?= 0

# Submodules
PWD = $(shell pwd)
SEQTK_ROOT = ${PWD}/src/htslib/

# Flags
CXX=g++
CXXFLAGS += -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${SEQTK_ROOT} -lz -lhts -Wl,-rpath,${SEQTK_ROOT}

# Additional flags for release/debug
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O9 -DNDEBUG
	#LDFLAGS += --static
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/htslib/*.c) $(wildcard src/htslib/htslib/*.h)
SVSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
TARGETS = .htslib src/svprops

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && cd ../../ && touch .htslib

src/svprops: .htslib $(SVSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/htslib && make clean
	rm -f $(TARGETS) $(TARGETS:=.o) .htslib
