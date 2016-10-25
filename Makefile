NAME := md

# directories
INC_DIRS := ./
BIN_DIR := ./bin

# flags
CFLAGS := -g -std=c11 -Wall
LDFLAGS := -lm

INCLUDE := $(addprefix -I, $(INC_DIRS))

C_FILES := $(wildcard *.c)
OBJ_C_FILES := ${C_FILES:%.c=%.o}

$(NAME): $(OBJ_C_FILES)
	mkdir -p $(BIN_DIR)
	$(CC) $^ $(LDFLAGS) -o $(BIN_DIR)/$@

%.o: %.c
	$(CC) $? $(CFLAGS) $(INCLUDE) -c -o $@

.PHONY: all clean 

clean:
	rm -rf *.o $(BIN_DIR)
