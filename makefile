PLATFORM=$(shell uname -s)

# Figure out bamtools, if we find the binary use it as the root or fallback to system libs
BAMTOOLS_BIN=$(shell command -v bamtools 2>/dev/null)
BAMTOOLS_PREFIX=$(patsubst %/bin/bamtools,%,$(BAMTOOLS_BIN))

# Guess include and lib dirs under that prefix
# don't use this if the user specified one
BAMTOOLS_INCLUDE ?= $(wildcard $(BAMTOOLS_PREFIX)/include/bamtools)
BAMTOOLS_LIB     ?= $(firstword $(wildcard $(BAMTOOLS_PREFIX)/lib64 $(BAMTOOLS_PREFIX)/lib))

GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iparallel-hashmap/parallel_hashmap/ -Icxxopts/include -Wno-unused-parameter `pkg-config --cflags zlib` -IMBG/src -Iconcurrentqueue

ifneq ($(BAMTOOLS_INCLUDE),)
CPPFLAGS+=-I$(BAMTOOLS_INCLUDE)
endif

ODIR=obj
BINDIR=bin
SRCDIR=src
LIBDIR=lib

LIBS=`pkg-config --libs zlib`
ifneq ($(BAMTOOLS_LIB),)
LIBS+=-L$(BAMTOOLS_LIB)
endif
LIBS+=-lbamtools

_DEPS = MatchIndex.h MinimizerIterator.h TwobitString.h ReadStorage.h UnitigKmerCorrector.h UnitigStorage.h ReadMatchposStorage.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = MatchIndex.o MinimizerIterator.o TwobitString.o ReadStorage.o UnitigKmerCorrector.o UnitigStorage.o ReadMatchposStorage.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))


ifeq ($(PLATFORM),Linux)
   LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++
else
   CPPFLAGS += -D_LIBCPP_DISABLE_AVAILABILITY
   LINKFLAGS = $(CPPFLAGS) $(LIBS) -lpthread -pthread -static-libstdc++
endif

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)
$(shell mkdir -p lib)

all: $(BINDIR)/matchchains $(BINDIR)/matchchains_batched $(BINDIR)/matchchains_index $(BINDIR)/matchchains_matchindex lib

lib: $(LIBDIR)/hifioverlapper.a

$(LIBDIR)/hifioverlapper.a: $(OBJ) $(DEPS)
	ar rvs $@ $(OBJ) $^

$(BINDIR)/matchchains: $(OBJ) $(ODIR)/matchchains.o MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/matchchains.o: $(SRCDIR)/matchchains.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(BINDIR)/matchchains_batched: $(OBJ) $(ODIR)/matchchains_batched.o MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/matchchains_batched.o: $(SRCDIR)/matchchains_batched.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(BINDIR)/matchchains_index: $(OBJ) $(ODIR)/matchchains_index.o MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/matchchains_index.o: $(SRCDIR)/matchchains_index.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(BINDIR)/matchchains_matchindex: $(OBJ) $(ODIR)/matchchains_matchindex.o MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/matchchains_matchindex.o: $(SRCDIR)/matchchains_matchindex.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

MBG/lib/mbg.a:
	$(MAKE) -C MBG lib

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	rm -f $(LIBDIR)/*
	$(MAKE) -C MBG clean
