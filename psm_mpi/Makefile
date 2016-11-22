CFLAGS:=-Wall -O1 -std=c99
LDLIBS:=-lm
PROJECT_NAME:=psm
SOURCE:=main.c
CC:=mpicc

.PHONY: clean

$(PROJECT_NAME): $(SOURCE)
	$(CC) $(CFLAGS) $(LDLIBS) -o $@ $<

clean:
	rm -rf *.o $(PROJECT_NAME)
