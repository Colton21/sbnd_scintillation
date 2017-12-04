CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)


run : PSD_perEvt
	@echo "Complied... PSD_perEvt"
PSD_perEvt : PSD_perEvt.o

	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^

