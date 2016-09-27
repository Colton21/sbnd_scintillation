CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)

run : libraryanalyze_light_histo
			@echo "Finished Compiling..."
			@echo "To run: ./libraryanalyze_light_histo"

libraryanalyze_light_histo : libraryanalyze_light_histo.o library_access.o utility_functions.o

	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^
