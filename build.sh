#/bin/sh

set -xe

CC=gcc
CFLAGS="-Wall -O2 -ggdb"
INC="-I./raylib/include"
LIBS="-lstdc++ -lm -L./raylib/lib -lraylib -ldl -lpthread"

#g++ $CFLAGS -c mtwist.c -o mtwist.o
DEPS="mtwist.o"

$CC $CFLAGS $INC sho-pimc.c $DEPS -o sho-pimc.exe $LIBS 
