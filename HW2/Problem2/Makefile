CC=g++
HEADER=-I$(CURDIR)/include
LIB=-lglut -lGLU -lGL
STD=c++11

foo: main.o RigidBody.o
	$(CC) -o foo main.o RigidBody.o $(LIB)

main.o: main.cpp RigidBody.h
	$(CC) -c main.cpp $(LIB) $(HEADER) -std=$(STD)

RigidBody.o: RigidBody.cpp RigidBody.h
	$(CC) -c RigidBody.cpp $(HEADER) -std=$(STD)

clean:
	rm -f *.o foo