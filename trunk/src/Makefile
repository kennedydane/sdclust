#################################################################
CC=g++
CFLAGS = -Wall -c
CFLAGS_DEBUG = -Wall -g -D debug_it_
CFLAGS_OPT = -O3 -Wall
CFLAGS_PROF = -Wall -pg -O3
OUTPUT = sd_clust

#################################################################
objects = clusters.o d1and2.o fasta.o globals.o match.o sd_clust.o
sources = clusters.cpp d1and2.cpp fasta.cpp globals.cpp match.cpp sd_clust.cpp
#################################################################
all: ../bin/sd_clust

../bin/sd_clust: $(objects)
	$(CC) $(objects) -o ../bin/$(OUTPUT)

clusters.o: clusters.h clusters.cpp
	$(CC) $(CFLAGS) clusters.cpp

d1and2.o: d1and2.h d1and2.cpp fasta.h globals.h
	$(CC) $(CFLAGS) d1and2.cpp

fasta.o: fasta.h fasta.cpp match.h globals.h d1and2.h
	$(CC) $(CFLAGS) fasta.cpp

globals.o: globals.h globals.cpp
	$(CC) $(CFLAGS) globals.cpp

match.o: match.h match.cpp clusters.h fasta.h globals.h d1and2.h
	$(CC) $(CFLAGS) match.cpp

sd_clust.o: sd_clust.cpp globals.h match.h fasta.h
	$(CC) $(CFLAGS) sd_clust.cpp

debug:
	$(CC) -g $(CFLAGS_DEBUG) $(sources) -o ../bin/$(OUTPUT).debug

opt:
	$(CC) $(CFLAGS_OPT) $(sources) -o ../bin/$(OUTPUT).opt

prof:
	$(CC) $(CFLAGS_PROF) $(sources) -o ../bin/$(OUTPUT).profile

clean:
	rm ../bin/$(OUTPUT)* $(objects)

