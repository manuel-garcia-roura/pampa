SRC_DIR = .
OBJ_DIR = .
LIB_DIR = ../../lib
BIN_DIR = ../../bin

SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES))

F90_FLAGS = -O3 -g -Wall -Werror
LD_FLAGS = -O3 -g -L$(LIB_DIR) -lpampa

LD_FLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LD_FLAGS += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

$(BIN_DIR)/pampa-f90: $(OBJ_FILES)
	mpif90 -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	mpif90 -c -o $@ $< $(F90_FLAGS)

clean:
	rm -f *.o
