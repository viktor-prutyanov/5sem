NAME := md_mpi

#compiler
OMPI_CXX := mpic++

# directories
INC_DIRS := ./
BIN_DIR := ./bin

# flags
OMPI_CPPFLAGS := -O2 -g -Wall
OMPI_LDFLAGS := -lm

INCLUDE := $(addprefix -I, $(INC_DIRS))

CPP_FILES := $(wildcard *.cpp)
OBJ_CPP_FILES := ${CPP_FILES:%.cpp=%.o}

$(NAME): $(OBJ_CPP_FILES)
	mkdir -p $(BIN_DIR)
	$(OMPI_CXX) $^ $(OMPI_LDFLAGS) -o $(BIN_DIR)/$@

%.o: %.cpp
	$(OMPI_CXX) $^ $(OMPI_CPPFLAGS) $(INCLUDE) -c -o $@

.PHONY: all clean 

clean:
	rm -rf *.o $(BIN_DIR)
