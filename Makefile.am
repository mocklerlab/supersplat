AUTOMAKE_OPTIONS = foreign subdir-objects

CXXFLAGS = -O3 -W -Wall -lrt -fopenmp -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0 $(SEQAN_CFLAGS)
LDFLAGS =

bin_PROGRAMS = supersplat
supersplat_SOURCES = src/main.cpp src/main.h src/Reads.cpp src/Reads.h src/Reference.cpp src/Reference.h src/Scheduler.cpp src/Scheduler.h src/Splat.cpp src/Splat.h src/Stack.cpp src/Stacker.cpp src/Stacker.h src/Stack.h src/StackParameters.cpp src/StackParameters.h src/stdafx.h src/SuperSplat.cpp
