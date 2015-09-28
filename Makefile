DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
SEQTK_ROOT = ${PWD}/src/htslib/

# Flags
CXX=g++
CXXFLAGS += -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${SEQTK_ROOT}

# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz
else
	LDFLAGS += -lhts -lz -Wl,-rpath,${SEQTK_ROOT}
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O9 -DNDEBUG
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SVSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
TARGETS = .htslib src/svprops src/sampleprops

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

src/svprops: .htslib $(SVSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/sampleprops: .htslib $(SVSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/htslib && make clean
	rm -f $(TARGETS) $(TARGETS:=.o) .htslib
