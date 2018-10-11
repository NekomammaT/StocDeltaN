MODEL = double_chaotic_conf
SOURCE = ../source
PHASECONF = conf

CXX := g++-8
CXXFLAGS := -std=c++11 -fopenmp -O3

ifeq ($(PHASECONF),conf)
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o

$(MODEL).o StocDeltaN_conf.o JacobiPDE_conf.o SRK32_conf.o: \
	$(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/StocDeltaN_conf.cpp \
	$(SOURCE)/JacobiPDE_conf.hpp $(SOURCE)/JacobiPDE_conf.cpp \
	$(SOURCE)/SRK32_conf.hpp $(SOURCE)/SRK32_conf.hpp
endif

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
