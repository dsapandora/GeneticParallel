#include "ParallelGenetics.h"

namespace ParallelGenetics{
	std::default_random_engine engine;
	std::uniform_int_distribution<int> parentDistribution;
	std::uniform_int_distribution<int>* crossoverDistributions;

	IndividualProperties::IndividualProperties(const std::vector<int>& genesPerChromosome){
		this->genesPerChromosome = genesPerChromosome;
	}

	std::vector<int> IndividualProperties::getProperties() const{
		return genesPerChromosome;
	}

	std::default_random_engine Individual::randomEngine = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

	Individual::Individual(const Individual& individual){
		this->chromosomes = individual.getChromosomes();
		this->evaluation = individual.getEvaluation();
		this->evaluated = individual.hasEvaluation();
	}

	Individual::Individual(const IndividualProperties& properties){
		std::vector<int> genesPerChromosome = properties.getProperties();
		std::uniform_int_distribution<int> distribution(0,1);
		for(int i=0; i<genesPerChromosome.size(); i++){
			chromosomes.push_back(std::vector<bool>());
			int numGenes = genesPerChromosome[i];
			for(int j=0; j<numGenes; j++){
				chromosomes[i].push_back(distribution(randomEngine));
			}
		}
		evaluated = false;
	}

	Individual::Individual(const std::vector<std::vector<bool>>& chromosomes){
		this->chromosomes = chromosomes;
		evaluated = false;
	}

	std::vector<std::vector<bool>> Individual::getChromosomes() const{
		return chromosomes;
	}

	std::vector<bool> Individual::getChromosomesAsString() const{
		std::vector<bool> genesString;

		for(int i=0; i<chromosomes.size(); i++){
			for(int j=0; j<chromosomes[i].size(); j++){
				genesString.push_back(chromosomes[i][j]);
			}
		}

		return genesString;
	}

	double Individual::getEvaluation() const{
		return evaluation;
	}

	void Individual::setEvaluation(double evaluation){
		this->evaluation = evaluation;
		evaluated = true;
	}

	bool Individual::hasEvaluation() const{
		return evaluated;
	}

	Population::Population(int numIndividuals, const IndividualProperties& properties){
		for(int i=0; i<numIndividuals; i++){
			population.push_back(Individual(properties));
		}
	}

	Population::Population(const std::vector<Individual>& individuals){
		population = individuals;
	}

	std::vector<Individual> Population::getPopulation() const{
		return population;
	}

	void Population::setPopulation(const std::vector<Individual>& population){
		this->population = population;
	}

	int Population::getNumIndividuals() const{
		return population.size();
	}

	Individual Population::getIndividual(int i) const{
		return population[i];
	}

	void Population::setIndividual(int i, const Individual& individual){
		population[i] = individual;
	}

	void Population::addIndividual(const Individual& individual){
		population.push_back(individual);
	}

	unsigned long long binToDec(std::vector<bool> vec){
		long long exp = 1;
		long long acum = 0;
		for(int i=vec.size()-1; i>=0; i--){
			if(vec[i]) acum += exp;
			exp *= 2;
		}

		return acum;
	}

	//Imprime a los individuos de la población
	void printPopulation(const Population& pop){
		std::vector<Individual> population = pop.getPopulation();

		for(int ind=0; ind<population.size(); ind++){
			std::vector<std::vector<bool>> chromosomes = population[ind].getChromosomes();
			std::cout << "\t";
			std::cout << ind+1 << ")\t";
			for(int i=0; i<chromosomes.size(); i++){
				for(int j=0; j<chromosomes[i].size(); j++){
					if(chromosomes[i][j]) std::cout << 1;
					else std::cout << 0;
				}
				std::cout << " ";
			}
			std::cout << std::endl;
		}
	}

	//Imprime datos sobre la generación actual
	void printData(int currentGeneration, const Population& pop, const Individual* best, int verbosity){
		std::cout << "Generation " << currentGeneration << ":" << std::endl;
		if(verbosity == PG_OUTPUT_COMPLETE){
			printPopulation(pop);
		}
		std::cout << std::endl;
		if(best != NULL){
			std::cout << "\tBest individual found so far: ";
			std::vector<std::vector<bool>> chromosomes = best->getChromosomes();
			for(int i=0; i<chromosomes.size(); i++){
				for(int j=0; j<chromosomes[i].size(); j++){
					if(chromosomes[i][j]) std::cout << 1;
					else std::cout << 0;
				}
				std::cout << " ";
			}
		} else {
			std::cout << "No best individual yet";
		}
		std::cout << std::endl << std::endl;
	}

	Population Select(const Population& pop, double (*evaluator)(const Individual&), Individual** best){
		std::vector<Individual> population = pop.getPopulation();		
		if(*best == NULL){
			*best = new Individual(population[0]);
		}
		for(int i=0; i<population.size(); i++){
			if(!population[i].hasEvaluation()){
				population[i].setEvaluation((*evaluator)(population[i]));
				if(population[i].getEvaluation() > (*best)->getEvaluation()){
					delete *best;
					*best = new Individual(population[i]);
				}
			}
		}
		shuffle(population.begin(), population.end(), engine);

		std::vector<Individual> newPopulation;
		if(population[0].getEvaluation() > population[population.size()-1].getEvaluation()){
			newPopulation.push_back(population[0]);
		} else {
			newPopulation.push_back(population[population.size()-1]);
		}
		for(int i=1; i<population.size(); i++){
			if( population[i-1].getEvaluation() > population[i].getEvaluation() ){
				newPopulation.push_back(population[i-1]);
			} else {
				newPopulation.push_back(population[i]);
			}
		}

		return Population(newPopulation);
	}

	Population Crossover(double crossoverProbability, const Population& pop){
		std::vector<Individual> oldPop = pop.getPopulation();
		std::vector<Individual> newPop;

		for(int i=0; i<oldPop.size()/2; i++){
			int first = parentDistribution(engine);
			int second = parentDistribution(engine);
			std::vector<std::vector<bool>> firstParent = oldPop[first].getChromosomes();
			std::vector<std::vector<bool>> secondParent = oldPop[second].getChromosomes();
			if( ((double)engine())/((double)engine.max()) <= crossoverProbability ){
				std::vector<std::vector<bool>> firstChild;
				std::vector<std::vector<bool>> secondChild;
				for(int c=0; c<firstParent.size(); c++){
					firstChild.push_back(std::vector<bool>());
					secondChild.push_back(std::vector<bool>());
					int crossoverPoint = crossoverDistributions[c](engine);
					for(int g=0; g<firstParent[c].size(); g++){
						if(g < crossoverPoint){
							firstChild[c].push_back(firstParent[c][g]);
							secondChild[c].push_back(secondParent[c][g]);
						} else{
							firstChild[c].push_back(secondParent[c][g]);
							secondChild[c].push_back(firstParent[c][g]);
						}
					}
				}
				newPop.push_back(Individual(firstChild));
				newPop.push_back(Individual(secondChild));
			} else{
				newPop.push_back(oldPop[first]);
				newPop.push_back(oldPop[second]);
			}
		}

		return Population(newPop);
	}

	Population Mutate(double mutationProbability, const Population& pop){
		std::vector<Individual> oldPop = pop.getPopulation();
		std::vector<Individual> newPop;
		bool wasMutated = false;

		for(int i=0; i<oldPop.size(); i++){
			std::vector<std::vector<bool>> chromosomes = oldPop[i].getChromosomes();
			for(int c=0; c<chromosomes.size(); c++){
				for(int g=0; g<chromosomes[c].size(); g++){
					if( ((double)engine())/((double)engine.max()) <= mutationProbability ){
						chromosomes[c][g] = !chromosomes[c][g];
						wasMutated = true;
					}
				}
			}
			if(wasMutated){
			newPop.push_back(Individual(chromosomes));
			} else {
				newPop.push_back(oldPop[i]);
			}
		}

		return Population(newPop);
	}

	Individual GeneticAlgorithmRun(int numGenerations, int numIndividuals, IndividualProperties individualProperties, double crossoverProbability, double mutationProbability, double (*evaluator)(const Individual&), int verbosity, Population** out){
		assert(numGenerations > 0);
		assert((mutationProbability>=0) && (mutationProbability<=1));
		assert((numIndividuals%2 == 0) && (numIndividuals >= 2));

		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		engine.seed(seed);
		parentDistribution = std::uniform_int_distribution<int>(0,numIndividuals-1);
		std::vector<int> genesPerChromosome = individualProperties.getProperties();
		crossoverDistributions = new std::uniform_int_distribution<int>[genesPerChromosome.size()];
		for(int i=0; i<genesPerChromosome.size(); i++){
			crossoverDistributions[i] = std::uniform_int_distribution<int>(1,genesPerChromosome[i]-1);
		}

		Population pop = Population(numIndividuals, individualProperties);
		Individual* bestSoFar = NULL;

		for(int g=1; g<=numGenerations; g++){
			pop = Select(pop, evaluator, &bestSoFar);
			pop = Crossover(crossoverProbability, pop);
			pop = Mutate(mutationProbability, pop);
			
			if(verbosity != PG_OUTPUT_NONE) printData(g,pop,bestSoFar,verbosity);
		}

		if(out != NULL) *out = new Population(pop);
		Individual best = Individual(bestSoFar->getChromosomes());
		delete [] crossoverDistributions;
		delete bestSoFar;
		return best;
	}
}