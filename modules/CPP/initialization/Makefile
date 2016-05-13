

objects = initialantigen_rmsd_main.o InitialAntigen.o
initialantigen_rmsd.out : $(objects)
	#g++ -L/usr/global/boost/1_51_0/lib -o initialantigen.out $(objects)
	g++ -o initialantigen_rmsd.out $(objects)
initialantigen_rmsd_main.o : initialantigen_rmsd_main.cpp InitialAntigen.h MOLECULES.h
	g++ -I/usr/global/boost/1_51_0/include -c initialantigen_rmsd_main.cpp -o initialantigen_rmsd_main.o -std=c++0x
InitialAntigen.o : InitialAntigen.cpp InitialAntigen.h MOLECULES.h ENERGY_FUNCTIONS.h DOCKING_FUNCTIONS.h
	g++  -I/usr/global/boost/1_51_0/include -c InitialAntigen.cpp -o InitialAntigen.o -std=c++0x

clean :
	rm -f $(objects)
