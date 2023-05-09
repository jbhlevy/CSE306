OBJS  = main.o
SOURCE  = main.cpp
HEADER  = 
OUT  = main
CC   = g++
CXXFLAGS = -g -O3 -std=c++11 -c -Wall -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
LFLAGS   =  -L/opt/homebrew/opt/llvm/lib/ -L/opt/homebrew/opt/libomp/lib -Xpreprocessor -fopenmp -lomp -O3

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.cpp
	$(CC) $(CXXFLAGS) main.cpp

clean:
	rm -f $(OBJS) $(OUT)
