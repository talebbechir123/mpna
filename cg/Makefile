CXX = g++
CXXFLAGS = -std=c++17 -O2

OBJ = CSRMatrix.o cg.o iterative_methods.o Laplacian2D.o main.o
EXE = solver

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE)

CSRMatrix.o: CSRMatrix.cpp CSRMatrix.h
	$(CXX) $(CXXFLAGS) -c CSRMatrix.cpp

cg.o: cg.cpp CSRMatrix.h
	$(CXX) $(CXXFLAGS) -c cg.cpp

iterative_methods.o: iterative_methods.cpp iterative_methods.h CSRMatrix.h
	$(CXX) $(CXXFLAGS) -c iterative_methods.cpp

Laplacian2D.o: Laplacian2D.cpp Laplacian2D.h CSRMatrix.h
	$(CXX) $(CXXFLAGS) -c Laplacian2D.cpp

main.o: main.cpp CSRMatrix.h iterative_methods.h Laplacian2D.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJ) $(EXE)
