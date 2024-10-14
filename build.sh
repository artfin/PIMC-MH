#/bin/sh

set -xe

CC=gcc
MPICC=mpicc
CFLAGS="-Wall -Wextra -O2 -ggdb"

INC="-I./raylib/include -I./raygui-4.0/ -I./raygui-4.0/src"

#INC="-I./raylib/include"
CLIENT_LIBS="-lstdc++ -lm"
SERVER_LIBS="-lstdc++ -lm -L./raylib/lib -l:libraylib.a -ldl -lpthread -lgsl -lgslcblas"
DEPS="mtwist.o"

#$CC $CFLAGS -c mtwist.c -o mtwist.o

#$CC $CFLAGS dummy_server.c -o ./dummy_server.exe 
$CC $CFLAGS $INC server.c $DEPS -o ./server.exe $SERVER_LIBS
$MPICC -DNO_VISUALIZE $CFLAGS $INC 1d_sho.c $DEPS -o 1d_sho.exe $CLIENT_LIBS 

#$CC $CFLAGS $INC m0_hear.c $DEPS -o m0_hear.exe $LIBS
#$CC $CFLAGS $INC 3d_sho.c $DEPS -o 3d_sho.exe $LIBS 

#$CC $CFLAGS $INC $INC_GUI animation_curve.c -o animate_curve.exe -L./raylib-5.0_linux_amd64/lib -l:libraylib.a -lm 
