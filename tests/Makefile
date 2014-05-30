CXX = g++ -std=gnu++11

MACOS = macos
ifeq ($(MACOS), macos)
	EXTRA_LIB_FLAGS = -framework Accelerate
endif

LIBCHEB = ~/Develop/cheb++/build/lib/libcheb.a

LIB_FLAGS = $(LIBCHEB) -larmadillo -lmx -lmex -lmat $(EXTRA_LIB_FLAGS)

OPT = -O2

#ARMAFINAL = -DARMA_NO_DEBUG

#WARN = -Wall
## Uncomment the above line to enable all compilation warings.

#BENCHMARK = -pg
## Uncomment the above line to enable bechmark mode.

HEADERS = -I../include
CXXFLAGS = $(BENCHMARK) $(WARN) $(ARMAFINAL) $(OPT) $(HEADERS)

all: test_D test_Boundary test_quadrature test_interpolation

test_D: test_D.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIB_FLAGS)

test_Boundary: test_Boundary.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIB_FLAGS)

test_quadrature: test_quadrature.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIB_FLAGS)

test_interpolation: test_interpolation.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIB_FLAGS)

.PHONY: clean
clean:
	rm -f test_D test_Boundary test_quadrature test_interpolation
