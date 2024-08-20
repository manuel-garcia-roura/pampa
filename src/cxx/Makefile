SRC_DIR = .
OBJ_DIR = .
LIB_DIR = ../../lib
BIN_DIR = ../../bin

SRC_FILES = $(wildcard $(SRC_DIR)/*.cxx)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))

CXX_FLAGS = -O3 -g -Wall -Werror -I..
LD_FLAGS = -O3 -g -L$(LIB_DIR) -lpampa

CXX_FLAGS += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
CXX_FLAGS += -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include
LD_FLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LD_FLAGS += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

$(BIN_DIR)/pampa-cxx: $(OBJ_FILES)
	mpic++ -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	mpic++ -c -o $@ $< $(CXX_FLAGS)

clean:
	rm -f *.o
