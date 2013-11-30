all:
	g++ -o GAfold GAfold.cpp GAevolve.cpp vienna/lib.a -O3