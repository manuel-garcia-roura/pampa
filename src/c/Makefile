SRC_DIR = .
OBJ_DIR = .
LIB_DIR = ../../lib
BIN_DIR = ../../bin

SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES))

C_FLAGS = -O3 -g -Wall -Werror -I..
LD_FLAGS = -O3 -g -L$(LIB_DIR) -lpampa

LD_FLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LD_FLAGS += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

$(BIN_DIR)/pampa-c: $(OBJ_FILES)
	mpic++ -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mpicc -c -o $@ $< $(C_FLAGS)

clean:
	rm -f *.o
