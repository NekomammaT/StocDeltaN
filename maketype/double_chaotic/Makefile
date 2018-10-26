MODEL = double_chaotic_conf
SOURCE = ../source
PHASECONF = conf

CXX := g++-8
CXXFLAGS := -std=c++11 -O3 -fopenmp

ifeq ($(PHASECONF),conf)
$(MODEL): $(MODEL).o $(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o \
	$(SOURCE)/StocDeltaN_conf.o $(SOURCE)/JacobiPDE_conf.o $(SOURCE)/SRK32_conf.o

$(MODEL).o: $(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/JacobiPDE_conf.hpp $(SOURCE)/SRK32_conf.hpp
StocDeltaN_conf.o: $(SOURCE)/StocDeltaN_conf.hpp $(SOURCE)/JacobiPDE_conf.hpp \
	$(SOURCE)/SRK32_conf.hpp $(SOURCE)/matplotlibcpp.hpp
JacobiPDE_conf.o: $(SOURCE)/JacobiPDE_conf.hpp
SRK32_conf.o: $(SOURCE)/SRK32_conf.hpp
endif

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
