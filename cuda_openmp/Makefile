CC = nvcc
FLAGS = -Xcompiler -fopenmp

default: cudm

cudm: cudm.o 
	$(CC) $(FLAGS) -o cudm cudm.o 

cudm.o: cudm.cu FeatureExtractor.h Sketch.h SketchIO.h
	$(CC) $(FLAGS) -c cudm.cu

clean:
	$(RM) cudm *.o *~
