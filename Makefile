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

# LIB_FLAGS = -larmadillo $(EXTRA_LIB_FLAGS)
LIB_FLAGS = $(EXTRA_LIB_FLAGS)
## On MacOS, we can use Armadillo without installation of Armadillo by using -framework Accelerate.
## On other platform, we should also add -DARMA_DONT_USE_WRAPPER -lblas -llapack

OPT = -O2

#ARMAFINAL = -DARMA_NO_DEBUG

#WARN = -Wall
## Uncomment the above line to enable all compilation warings.

#BENCHMARK = -pg
## Uncomment the above line to enable bechmark mode.

LIB_HEADERS = -Iinclude -I/Users/lyx/SynologyDrive/Develop/others/armadillo-9.900.2/include

CXXFLAGS = $(BENCHMARK) $(WARN) $(ARMAFINAL) $(OPT) $(LIB_HEADERS)

MKDIR_P = mkdir -p

_OBJS = Cheb.o \
		Boundary.o

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
