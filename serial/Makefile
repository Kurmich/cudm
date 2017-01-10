CC = g++
FLAGS = -fopenmp

default: cudm

cudm: cudm.o 
	$(CC) -o cudm cudm.o 

cudm.o: cudm.cpp FeatureExtractor.h Sketch.h SketchIO.h
	$(CC) $(FLAGS) -c cudm.cpp

clean:
	$(RM) cudm *.o *~
