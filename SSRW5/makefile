.PHONY: SSRW5 SSRW5.o
SSRW5: SSRW5.o SFMT.o
		g++ -g -std=c++11  -pthread -o SSRW5 SSRW5.o SFMT.o -ligraph
SSRW5.o: SFMT/SFMT.h 
		g++ -g -std=c++11  -pthread -c SSRW5.cpp
SFMT.o: SFMT/SFMT.h
	gcc  -c SFMT/SFMT.c
clean:
		-rm  *.o
