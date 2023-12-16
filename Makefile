HybridGen40:
	g++ -O3 -std=c++11 -fopenmp mpiwithopenmp.cpp -o mwo40 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

HybridGen80:
	g++ -O3 -std=c++11 -fopenmp mpiwithopenmp.cpp -o mwo0 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

HybridGen160:
	g++ -O3 -std=c++11 -fopenmp mpiwithopenmp.cpp -o mwo160 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

Sub40:
	./Sub40.sh

Sub80:
	./Sub80.sh

Sub160:
	./Sub160.shu