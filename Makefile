DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -fno-strict-aliasing
LDFLAGS += -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib

# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG -D__STDC_LIMIT_MACROS
endif
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/svprops src/sampleprops
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

src/svprops: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/sampleprops: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	cd src/htslib && make clean
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}
