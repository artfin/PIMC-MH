#/bin/sh

set -xe

CC=gcc
MPICC=mpicc
CFLAGS="-Wall -Wextra -O2 -ggdb"
INC="-I./raylib/include"
LIBS="-lstdc++ -lm -L./raylib/lib -lraylib -ldl -lpthread -lgsl -lgslcblas"
DEPS="mtwist.o"

$CC $CFLAGS -c mtwist.c -o mtwist.o

$CC $CFLAGS $INC server.c $DEPS -o ./server.exe $LIBS
$MPICC $CFLAGS $INC 1d_sho.c $DEPS -o 1d_sho.exe $LIBS 

#$CC $CFLAGS $INC m0_hear.c $DEPS -o m0_hear.exe $LIBS
#$CC $CFLAGS $INC 3d_sho.c $DEPS -o 3d_sho.exe $LIBS 
