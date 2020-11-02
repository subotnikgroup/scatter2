CXX=g++
FLAG=-O3 -g
THREADFLAG=-fopenmp -lpthread
EXFLAG=

all: scatter

scatter: scatter.cpp potential.h angularpot.hpp linearpot.hpp printmat.hpp
	python3 preprocfac.py angularpot.hpp linearpot.hpp scatter.cpp > scatter1.cpp
	$(CXX) -o scatter scatter1.cpp $(FLAG) $(THREADFLAG) $(EXFLAG) -larmadillo
	rm scatter1.cpp

clean:
	rm scatter

