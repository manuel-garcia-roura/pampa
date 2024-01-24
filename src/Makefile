all:
	rm -rf ../bin
	mkdir ../bin
	mpic++ -g *.cxx -o ../bin/pampa \
		-I/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include -L/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib -lpetsc \
		-I/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/include -L/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib -lslepc
