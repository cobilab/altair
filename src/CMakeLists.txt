cmake_minimum_required(VERSION 2.8.4)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Ofast -march=native -flto -ffast-math -fomit-frame-pointer -funroll-loops -lpthread -lm -w " )
IF(UNIX)
 link_libraries(m)
ENDIF(UNIX)

project (ALTAIR) 
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_BUILD_TYPE Release)

add_executable (AltaiR altair.c rw.c ncd.c nc.c fc.c af.c ma.c kmer.c stats.c dist.c cm.c chash.c carray.c rand.c tolerant.c pmodels.c levels.c cache.c math.c dna.c keys.c alphabet.c buffer.c param.c strings.c mem.c msg.c common.c threads.h defs.h) 

TARGET_LINK_LIBRARIES(AltaiR pthread)
