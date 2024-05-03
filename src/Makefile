SRC_DIR = .
OBJ_DIR = .
BIN_DIR = ../bin

PETSC_DIR = /lib/petscdir/petsc3.12
SLEPC_DIR = /lib/slepcdir/slepc3.12
PETSC_ARCH = x86_64-linux-gnu-real

EIGEN_DIR = /usr/include/eigen3
WITH_EIGEN = true

METIS_DIR = /usr/local
WITH_METIS = false

SRC_FILES = $(wildcard $(SRC_DIR)/*.cxx)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))

CXX_FLAGS = -O3 -Wall -Werror
LD_FLAGS = -O3

CXX_FLAGS += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
CXX_FLAGS += -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include
LD_FLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
LD_FLAGS += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

ifeq ($(WITH_EIGEN), true)
CXX_FLAGS += -I$(EIGEN_DIR) -DWITH_EIGEN
endif

ifeq ($(WITH_METIS), true)
CXX_FLAGS += -I$(METIS_DIR)/include -DWITH_METIS
LD_FLAGS += -L$(METIS_DIR)/lib -lmetis
endif

$(BIN_DIR)/pampa: $(OBJ_FILES)
	rm -rf $(BIN_DIR)
	mkdir $(BIN_DIR)
	mpic++ -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	mpic++ -c -o $@ $< $(CXX_FLAGS)

clean:
	rm -f *.o
