OBJS  = main.o lbfgs.o
SOURCE  = main.cpp lbfgs.cpp
HEADER  = 
OUT  = main
CC   = g++
CXXFLAGS = -g -O3 -std=c++11 -c -Wall -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
LFLAGS   =  -L/opt/homebrew/opt/llvm/lib/ -L/opt/homebrew/opt/libomp/lib -Xpreprocessor -fopenmp -lomp -O3

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.cpp
	$(CC) $(CXXFLAGS) main.cpp

lbfgs.o: lbfgs.cpp
	$(CC) $(CXXFLAGS) lbfgs.cpp

clean:
	rm -f $(OBJS) $(OUT)
	rm -f *.png
