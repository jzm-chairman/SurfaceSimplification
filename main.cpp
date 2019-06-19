#include "Simplifier.hpp"

int main(int argc, char* argv[])
{
	if (argc != 4) { printf("Usage: %s <input file> <output file> <simplification ratio>\n", argv[0]); exit(0); }
	clock_t start = clock();
	Simplifier* simp = new Simplifier(atof(argv[3]), 0);
	simp->read(argv[1]);
	simp->buildHeap();
	simp->runSimp();
	simp->write(argv[2]);
	clock_t finish = clock();
	printf("Model %s, Elapsed Time %fs\n", argv[1], (double)(double(finish - start) / CLOCKS_PER_SEC));
	return 0;
}