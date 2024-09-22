#/bin/sh

set -xe

CC=gcc
CFLAGS="-Wall -Wextra -O2 -ggdb"
INC="-I./raylib/include"
LIBS="-lstdc++ -lm -L./raylib/lib -lraylib -ldl -lpthread -lgsl -lgslcblas"
DEPS="mtwist.o"

gcc $CFLAGS -c mtwist.c -o mtwist.o

$CC $CFLAGS $INC 1d_sho.c $DEPS -o 1d_sho.exe $LIBS 
#$CC $CFLAGS $INC 3d_sho.c $DEPS -o 3d_sho.exe $LIBS 
