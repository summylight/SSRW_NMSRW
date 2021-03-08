.PHONY: NMSRW5 NMSRW5.o
NMSRW5: NMSRW5.o SFMT.o
		g++ -g -std=c++11  -pthread -o NMSRW5 NMSRW5.o SFMT.o -ligraph
NMSRW5.o: SFMT/SFMT.h 
		g++ -g -std=c++11  -pthread -c NMSRW5.cpp
SFMT.o: SFMT/SFMT.h
	gcc  -c SFMT/SFMT.c
clean:
		-rm  *.o
