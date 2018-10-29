CC=g++
HEADER=-I$(CURDIR)/include
LIB=-lglut -lGLU -lGL
STD=c++11

foo: main.o RigidBody.o Contact.o
	$(CC) -o foo main.o RigidBody.o Contact.o $(LIB)

main.o: main.cpp Contact.h
	$(CC) -c main.cpp $(LIB) $(HEADER) -std=$(STD)

RigidBody.o: RigidBody.cpp RigidBody.h
	$(CC) -c RigidBody.cpp $(HEADER) -std=$(STD)

Contact.o: Contact.cpp Contact.h RigidBody.h
	$(CC) -c Contact.cpp $(HEADER) -std=$(STD)

clean:
	rm -f *.o foo