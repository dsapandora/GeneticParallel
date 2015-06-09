#include <iomanip>
#include "ParallelGenetics.h"

using namespace std;
using namespace ParallelGenetics;

double evalFunction(const Individual& individual){
	vector<vector<bool> > chromosomes = individual.getChromosomes();
	static double intervalSize = ((double)10)/pow(2,32);
	unsigned long long decX = binToDec(chromosomes[0]);
	double xVal = -5.0 + intervalSize*decX;
	double f = exp(xVal*cos(xVal));

	return f;
}

double pozos(const Individual& individual){
	double costoPlataforma[] = {31.7500, 31.0454, 34.1459, 36.7451, 58.0218, 24.2361,31.5000, 32.1231, 32.1231, 42.9861, 61.4721, 39.7500,67.1053, 69.7658, 99.2318, 49.4890, 23.5000, 367.3722};
	int pozosCubiertosPorPlataforma[] = {7, 56, 448, 416, 1984, 1536, 98304, 6144, 24576,155648, 30720, 917504, 1015808, 99840, 511, 31,786432, 1048575};
	double costoTodasPlataformas = 0;

	int indDec = binToDec(individual.getChromosomes()[0]);
	double value = 0;
	int pozosCubiertos = 0;

	for(int i=0; i<18; i++){
		if(indDec & 1){
			value += costoPlataforma[i];
			pozosCubiertos |= pozosCubiertosPorPlataforma[i];
		}
		indDec >>= 1;
		costoTodasPlataformas += costoPlataforma[i];
	}

	int pozosSinPlataforma = 0;
	if(pozosCubiertos != 1048575){
		for(int i=0; i<20; i++){
			if( (pozosCubiertos >> i) & 1 ){
				pozosSinPlataforma += 1;
			}
		}
	}

	value += costoTodasPlataformas*pozosSinPlataforma;
	value = costoTodasPlataformas/value;

	if(pozosCubiertos == 0){
		value = -costoTodasPlataformas*100;
	}

	return value;
}


int main(){
	/* e^(x*cos(x)) */
	vector<int> genesPerChromosome;
	genesPerChromosome.push_back(32);
	IndividualProperties properties = IndividualProperties(genesPerChromosome);
	Individual best = GeneticAlgorithmRun(1000,2500,properties,0.7,0.0625,&evalFunction,PG_OUTPUT_NONE);

	unsigned long long decX = binToDec(best.getChromosomes()[0]);
	static double intervalSize = ((double)10)/pow(2,32);
	double xVal = -5.0 + intervalSize*decX;
	cout << fixed << setprecision(7) 
		<< "Mejor solucion encontrada: f(" << xVal << ") = " << evalFunction(best) << endl;

	
	/* Still not work, but someday will do :)
	vector<int> prop2;
	prop2.push_back(16);
	IndividualProperties properties = IndividualProperties(prop2);
	Individual best = runGeneticAlgorithm(100,1000,properties,0.7,0.055,&pozos,PG_OUTPUT_NONE);

	cout << "Plataformas usadas en la mejor solucion encontrada: ";
	unsigned long long num = binToDec(best.getChromosomes()[0]);
	for(int i=0; i<18; i++){
		if( (num >> i) & 1 ) cout << i+1 << " ";
	}*/

	return 0;
}