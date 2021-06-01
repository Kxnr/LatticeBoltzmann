CX=g++
CFLAGS=-O3 -Werror -Wpedantic -std=c++17 

ODIR=obj
IDIR=include
SDIR=src

IFLAGS=-Iinclude -Ilib
LFLAGS=-larmadillo

_OBJ=lbm.o
_DEPS=lbm.h
_SRC=main.cpp lbm.cpp

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
SRC = $(patsubst %,$(SDIR)/%,$(_SRC))


all: main 

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) $(CFLAGS) $(IFLAGS) -c -o $@ $<

main: $(OBJ)
	$(CXX) $(CFLAGS) $(IFLAGS) -o $@.out $^ $(LFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o 
