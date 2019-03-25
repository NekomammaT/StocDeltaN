MODEL = chaotic_conf
SOURCE = ../source

CXX := g++-8
CXXFLAGS := -std=c++11 -fopenmp

$(MODEL): $(MODEL).o $(SOURCE)/JacobiPDE.o 
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/JacobiPDE.o 

$(MODEL).o: $(SOURCE)/JacobiPDE.hpp 
JacobiPDE_conf.o: $(SOURCE)/JacobiPDE.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o	
	$(RM) $(MODEL)
