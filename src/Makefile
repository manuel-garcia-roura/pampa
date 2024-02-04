SRC_DIR := .
OBJ_DIR := .
BIN_DIR := ../bin

PETSC_DIR := /lib/petscdir/petsc3.12/x86_64-linux-gnu-real
SLEPC_DIR := /lib/slepcdir/slepc3.12/x86_64-linux-gnu-real

SRC_FILES := $(wildcard $(SRC_DIR)/*.cxx)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))

CXX_FLAGS := -O3 -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include
LD_FLAGS := -O3 -L$(PETSC_DIR)/lib -lpetsc -L$(SLEPC_DIR)/lib -lslepc

$(BIN_DIR)/pampa: $(OBJ_FILES)
	rm -rf $(BIN_DIR)
	mkdir $(BIN_DIR)
	mpic++ -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	mpic++ -c -o $@ $< $(CXX_FLAGS)

clean:
	rm -f *.o
