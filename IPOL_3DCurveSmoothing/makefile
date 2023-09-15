#   Copyright (c) 2010-2019
#   IMAGING TECHNOLOGY CENTER (CTIM - www.ctim.es) <lalvarez@dis.ulpgc.es>
#   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
#   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en


#Option to use OpenMP (use make command with OMP=1)
ifdef OMP
COPT	= -O3 -fopenmp
else
COPT	= -O3
endif

#Option to generate coverage test (use make command with COV=1)
ifdef COV
LCOV = --coverage
endif

#Compiler and linker options
CFLAGS	= -Wall -Wextra -std=c++98 $(COPT) $(LCOV)
LDFLAGS += -lm $(LCOV)

AUXILIARY_DIR = auxiliary/
MAIN_DIR      = main/

#Option to use the version of gcc installed with Hombrew in MacOS
#(modify the line if needed)
ifdef APPLE
CXX = /usr/local/bin/g++-9
endif

3DCurveSmoothing:\
		main3DCurveSmoothing.o curve3D.o obj3D.o
	$(CXX) $(CFLAGS) -o\
  3DCurveSmoothing main3DCurveSmoothing.o curve3D.o obj3D.o $(LDFLAGS)

main3DCurveSmoothing.o:\
	$(MAIN_DIR)main3DCurveSmoothing.cpp
	$(CXX) $(CFLAGS) -c $(MAIN_DIR)main3DCurveSmoothing.cpp

curve3D.o: $(MAIN_DIR)curve3D.cpp
	$(CXX) $(CFLAGS) -c $(MAIN_DIR)curve3D.cpp

obj3D.o: $(AUXILIARY_DIR)obj3D.cpp
	$(CXX) $(CFLAGS) -c $(AUXILIARY_DIR)obj3D.cpp

.PHONY: clean
clean:
	@echo "Cleaning files..."
	$(RM) -f *.o *.gcno *.gcda 3DCurveSmoothing
