MODEL = double_chaotic_conf
SOURCE = ../source

CXX := g++-8
CXXFLAGS := -std=c++11 -O3 -fopenmp

$(MODEL): $(MODEL).o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/JacobiPDE.o $(SOURCE)/SRK32.o

$(MODEL).o: $(SOURCE)/JacobiPDE.hpp $(SOURCE)/SRK32.hpp
JacobiPDE.o: $(SOURCE)/JacobiPDE.hpp
SRK32.o: $(SOURCE)/SRK32.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
