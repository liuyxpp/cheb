BASEDIR = build
LOUTDIR = $(BASEDIR)/lib
LOUT = $(LOUTDIR)/libcheb.a
ODIR = $(BASEDIR)/src
SDIR = src

CXX = g++ -std=gnu++11

MACOS = macos
ifeq ($(MACOS), macos)
	EXTRA_LIB_FLAGS = -framework Accelerate
endif

LIB_FLAGS = -larmadillo $(EXTRA_LIB_FLAGS)

OPT = -O2

#ARMAFINAL = -DARMA_NO_DEBUG

#WARN = -Wall
## Uncomment the above line to enable all compilation warings.

#BENCHMARK = -pg
## Uncomment the above line to enable bechmark mode.

LIB_HEADERS = -Iinclude

CXXFLAGS = $(BENCHMARK) $(WARN) $(ARMAFINAL) $(OPT) $(LIB_HEADERS)

MKDIR_P = mkdir -p

_OBJS = Cheb.o \
		Boundary.o \
		Etdrk4.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LOUT) : $(OBJS)
	ar -r $(LOUT) $^

.PHONY: all
all: $(LOUT)

.PHONY: lib
lib: $(LOUT)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o $(LOUT)
