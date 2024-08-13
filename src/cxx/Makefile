SRC_DIR = .
OBJ_DIR = .
LIB_DIR = ../../lib
BIN_DIR = ../../bin

EIGEN_DIR = /usr/include/eigen3
WITH_EIGEN = true

SRC_FILES = $(wildcard $(SRC_DIR)/*.cxx)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))

CXX_FLAGS = -O3 -g -Wall -Werror -I..
LD_FLAGS = -O3 -g -L$(LIB_DIR) -lpampa

CXX_FLAGS += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
CXX_FLAGS += -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include
LD_FLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LD_FLAGS += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

ifeq ($(WITH_EIGEN), true)
CXX_FLAGS += -I$(EIGEN_DIR) -DWITH_EIGEN
endif

$(BIN_DIR)/pampa: $(OBJ_FILES)
	rm -rf $(BIN_DIR)
	mkdir $(BIN_DIR)
	mpic++ -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	mpic++ -c -o $@ $< $(CXX_FLAGS)

clean:
	rm -f *.o
