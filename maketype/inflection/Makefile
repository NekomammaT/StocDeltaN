MODEL = inflection_logp
SOURCE = ../source
PHASECONF = phase

CXX := g++-8
CXXFLAGS := -std=c++11 -fopenmp -O2

ifeq ($(PHASECONF),conf)
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o

$(MODEL).o StocDeltaN_conf.o JacobiPDE_conf.o SRK32_conf.o: \
	$(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/StocDeltaN_conf.cpp \
	$(SOURCE)/JacobiPDE_conf.hpp $(SOURCE)/JacobiPDE_conf.cpp \
	$(SOURCE)/SRK32_conf.hpp $(SOURCE)/SRK32_conf.cpp
else
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN.o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN.o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o

$(MODEL).o StocDeltaN.o JacobiPDE.o SRK32.o: \
	$(SOURCE)/StocDeltaN.hpp $(SOURCE)/StocDeltaN.cpp \
	$(SOURCE)/JacobiPDE.hpp $(SOURCE)/JacobiPDE.cpp \
	$(SOURCE)/SRK32.hpp $(SOURCE)/SRK32.cpp
endif

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
