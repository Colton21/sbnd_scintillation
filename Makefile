CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)

run : libraryanalyze_light_histo
	./libraryanalyze_light_histo

libraryanalyze_light_histo : libraryanalyze_light_histo.o library_access.o utility_functions.o

#g++ -o libraryanalyze_light_histo LOADLIBRARYFROMFILE.o libraryanalyze_light_histo.o ${LIBS}
	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^

#LOADLIBRARYFROMFILE.o : LOADLIBRARYFROMFILE.cc LOADLIBRARYFROMFILE.h
#	g++ ${CXXFLAGS} -o LOADLIBRARYFROMFILE.o -c LOADLIBRARYFROMFILE.cc

#libraryanalyze_light_histo.o : libraryanalyze_light_histo.cc LOADLIBRARYFROMFILE.h
#	g++ ${CXXFLAGS} -o libraryanalyze_light_histo.o -c libraryanalyze_light_histo.cc
