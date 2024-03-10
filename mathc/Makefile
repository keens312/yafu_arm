CC =gcc 
CFLAGS= -g  -march=armv8-a -std=gnu99 -fPIC
SRC = mathc.c
OBJS =$(SRC:.c=$(OBJ_EXT))
HEAD= mathc.h 

all: $(OBJS)
	rm -f  libmathc.a
	ar r libmathc.a $(OBJS)
	ranlib libmathc.a
clean:
	rm -f $(OBJS)
%$(OBJ_EXT):%.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<
