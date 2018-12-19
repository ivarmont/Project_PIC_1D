OBJS = 1D_PIC.o Particle.o Domain.o Functions.o
CC = g++ 
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
LIBS = 

main :	$(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o main

1D_PIC.o : 1D_PIC.cpp Particle.h 
	$(CC) $(CFLAGS) 1D_PIC.cpp $(LIBS)

Domain.o : Domain.cpp Domain.h Functions.h
	$(CC) $(CFLAGS)  Domain.cpp $(LIBS)

Functions.o : Functions.cpp Functions.h
	$(CC) $(CFLAGS)  Functions.cpp $(LIBS)

Particle.o : Particle.cpp Particle.h 
	$(CC) $(CFLAGS)  Particle.cpp $(LIBS) 




clean:
	\rm *.o *~ main

