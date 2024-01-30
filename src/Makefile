SRC_DIR := .
OBJ_DIR := .
SRC_FILES := $(wildcard $(SRC_DIR)/*.cxx)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDFLAGS := -O3 -L/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib -lpetsc -L/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib -lslepc
CXXFLAGS := -O3 -I/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include -I/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/include

../bin/pampa: $(OBJ_FILES)
	rm -rf ../bin
	mkdir ../bin
	mpic++ -g -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	mpic++ -g -c -o $@ $< $(CXXFLAGS)

clean:
	rm *.o
