CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)


run : PSD_perEvt_convolution
	@echo "Complied... PSD_perEvt_convolution"
PSD_perEvt_convolution : PSD_perEvt_convolution.o

	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^

