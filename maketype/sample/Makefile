MODEL = double_chaotic_conf
SOURCE = ../source
PHASECONF = conf

CXX := g++
CXXFLAGS := -std=c++11 -O3 #-fopenmp

ifeq ($(PHASECONF),conf)
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o

$(MODEL).o: $(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/JacobiPDE_conf.hpp $(SOURCE)/SRK32_conf.hpp
StocDeltaN_conf.o: $(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/JacobiPDE_conf.hpp \
	$(SOURCE)/SRK32_conf.hpp $(SOURCE)/matplotlibcpp.hpp
JacobiPDE_conf.o: $(SOURCE)/JacobiPDE_conf.hpp
SRK32_conf.o: $(SOURCE)/SRK32_conf.hpp
else
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN.o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN.o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o

$(MODEL).o: $(SOURCE)/StocDeltaN.hpp $(SOURCE)/JacobiPDE.hpp $(SOURCE)/SRK32.hpp
StocDeltaN_conf.o: $(SOURCE)/StocDeltaN.hpp $(SOURCE)/JacobiPDE.hpp \
	$(SOURCE)/SRK32.hpp $(SOURCE)/matplotlibcpp.hpp
JacobiPDE_conf.o: $(SOURCE)/JacobiPDE.hpp
SRK32_conf.o: $(SOURCE)/SRK32.hpp
endif

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
