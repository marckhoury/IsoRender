TARGET = IsoRender
CC = g++
OBJECTS = ijkmcube_datastruct.o ijkmcubeIO.o ijkmcube_sub.o ijkoctree.o ijktable.o ijkmcube_extract.o ijkmcube.o ijkmcube_util.o ijksnapmc.o ijkxitIO.o Vector.o main.o
LIBS = -lGL -lglut -lGLU -lglui -lz -lexpat libITKNrrdIO.a libz.so

IsoRender: $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)

main.o: main.cpp
	$(CC) -c main.cpp

Vector.o: Vector.cpp Vector.h
	$(CC) -c Vector.cpp

ijkmcube_datastruct.o: ijkmcube_datastruct.cxx ijkmcube_datastruct.h
	$(CC) -c ijkmcube_datastruct.cxx

ijkmcubeIO.o: ijkmcubeIO.cxx ijkmcubeIO.h
	$(CC) -c ijkmcubeIO.cxx

ijkmcube_sub.o: ijkmcube_sub.cxx ijkmcube_sub.h
	$(CC) -c ijkmcube_sub.cxx

ijkoctree.o: ijkoctree.cxx ijkoctree.h
	$(CC) -c ijkoctree.cxx

ijktable.o: ijktable.cxx ijktable.h
	$(CC) -c ijktable.cxx

ijkmcube_extract.o:  ijkmcube_extract.cxx  ijkmcube_extract.h
	$(CC) -c ijkmcube_extract.cxx

ijkmcube.o: ijkmcube.cxx ijkmcube.h
	$(CC) -c ijkmcube.cxx

ijkmcube_util.o: ijkmcube_util.cxx ijkmcube_util.h
	$(CC) -c ijkmcube_util.cxx

ijksnapmc.o: ijksnapmc.cxx ijksnapmc.h
	$(CC) -c ijksnapmc.cxx

ijkxitIO.o: ijkxitIO.cxx ijkxitIO.h
	$(CC) -c ijkxitIO.cxx

clean:
	rm *.o IsoRender *~
