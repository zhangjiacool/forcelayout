CFLAGS=-I. -std=gnu99 -O3 -march=native -fno-math-errno -fno-trapping-math -fno-signed-zeros -ffinite-math-only
#CFLAGS=-I. -std=gnu99 -g -march=native
LDFLAGS=-ljansson -lm -lpthread
OBJS=world.o force.o adjust.o sparsify.o worker.o

all: forcelayout

forcelayout: $(OBJS)
	gcc $(LDFLAGS) -o forcelayout $(OBJS)

world.o: world.c world.h force.h adjust.h sparsify.h
	gcc -c $(CFLAGS) world.c

force.o: force.c force.h world.h adjust.h worker.h
	gcc -c $(CFLAGS) force.c

adjust.o: adjust.c adjust.h world.h worker.h
	gcc -c $(CFLAGS) adjust.c

sparsify.o: sparsify.c world.h worker.h
	gcc -c $(CFLAGS) sparsify.c

worker.o: worker.c worker.h
	gcc -c $(CFLAGS) worker.c

clean:
	rm -f $(OBJS) forcelayout
