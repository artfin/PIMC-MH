#/bin/sh

set -xe

CC=gcc
MPICC=mpicc
MPICXX=mpic++
CFLAGS="-Wall -Wextra -O2 -ggdb"

INC="-I./raylib/include -I./raygui/ -isystem./raygui/src"
INC_HEP="-I/home/artfin/Desktop/lib/hep-mc-0.7/include/"
#INC="-I./raylib/include"
CLIENT_LIBS="-lstdc++ -lm -L./raylib/lib -l:libraylib.a"
SERVER_LIBS="-lstdc++ -lm -L./raylib/lib -l:libraylib.a -ldl -lpthread -lgsl -lgslcblas"
DEPS="mtwist.o"

$CC $CFLAGS -c mtwist.c -o mtwist.o

$CC $CFLAGS ./tools/kahan.c -o ./tools/kahan.exe 
$CC $CFLAGS ./tools/dummy_server.c -o ./tools/dummy_server.exe 

$CC $CFLAGS $INC server.c $DEPS -o ./server.exe $SERVER_LIBS
$MPICC $CFLAGS $INC 1d_sho.c $DEPS -o 1d_sho.exe $CLIENT_LIBS 

$MPICC $CFLAGS $INC m0_hear.c $DEPS -o m0_hear.exe $CLIENT_LIBS
$CC $CFLAGS $INC 3d_sho.c $DEPS -o 3d_sho.exe $LIBS 

$MPICXX $CFLAGS $INC_HEP m0_hep.cpp -o m0_hep.exe

#$CC test.c binn_test.c -o ./test.exe
#$CC $CFLAGS $INC $INC_GUI animation_curve.c -o animate_curve.exe -L./raylib-5.0_linux_amd64/lib -l:libraylib.a -lm 
