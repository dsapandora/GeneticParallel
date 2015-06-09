#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cassert>
#ifdef max(a,b)
	#undef max(a,b)
#endif

namespace ParallelGenetics{
	enum verbosity {PG_OUTPUT_NONE, PG_OUTPUT_BASIC, PG_OUTPUT_COMPLETE};

	class IndividualProperties{
	private:
		std::vector<int> genesPerChromosome;
	public:
		IndividualProperties() {};
		IndividualProperties(const std::vector<int>& genesPerChromosome);
		std::vector<int> getProperties() const;
	};

	class Individual{
	private:
		std::vector<std::vector<bool> > chromosomes;
		double evaluation;
		bool evaluated;
		static std::default_random_engine randomEngine;
	public:
		Individual() {};
		Individual(const Individual& individual);
		Individual(const IndividualProperties& properties);
		Individual(const std::vector<std::vector<bool> >& chromosomes);
		std::vector<std::vector<bool> > getChromosomes() const;
		std::vector<bool> getChromosomesAsString() const;
		double getEvaluation() const;
		void setEvaluation(double evaluation);
		bool hasEvaluation() const;
	};

	class Population{
	private:
		std::vector<Individual> population;
	public:
		Population() {};
		Population(int numIndividuals, const IndividualProperties& properties);
		Population(const std::vector<Individual>& individuals);
		std::vector<Individual> getPopulation() const;
		void setPopulation(const std::vector<Individual>& population);
		int getNumIndividuals() const;
		Individual getIndividual(int i) const;
		void setIndividual(int i, const Individual& individual);
		void addIndividual(const Individual& individual);
	};

	unsigned long long binToDec(std::vector<bool>);

	Individual GeneticAlgorithmRun(int numGenerations,
		int numIndividuals,
		IndividualProperties individualProperties,
		double crossoverProbability,
		double mutationProbability,
		double (*evaluator)(const Individual&),
		int verbosity=0,
		Population** out=NULL);
}