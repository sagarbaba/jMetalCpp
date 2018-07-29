//Required Headers
#include <Problem.h>
#include <Algorithm.h>
#include <Solution.h>
#include <Operator.h>
#include <ProblemFactory.h>
#include <MOALCPSO.h>
#include <PolynomialMutation.h>
#include <iostream>
#include <time.h>

int main(int argc, char ** argv) {

  clock_t t_ini, t_fin;

  Problem   * problem;   // The problem to solve
  Algorithm * algorithm; // The algorithm to use
  Operator  * mutation;  // "Turbulence" operator

  map<string, void *> parameters; // Operator parameters


  //TODO: QualityIndicator * indicators; // Object to get quality indicators

  if (argc>=2) {
    problem = ProblemFactory::getProblem(argc, argv);
    cout << "Selected problem: " << problem->getName() << endl;
  } else {
    cout << "No problem selected." << endl;
    cout << "Default problem will be used: Kursawe" << endl;
    //cout<<"jhuvv";
    problem = ProblemFactory::getProblem(const_cast<char *>("Kursawe"));
    //cout<<"jhu";
  }

  algorithm = new MOALCPSO(problem);
  // Algorithm parameters
  int swarmSizeValue = 100;
  int archiveSizeValue = 100;
  int maxIterationsValue = 250;
  algorithm->setInputParameter("swarmSize",&swarmSizeValue);
  algorithm->setInputParameter("archiveSize",&archiveSizeValue);
  algorithm->setInputParameter("maxIterations",&maxIterationsValue);

  // Mutation operator
  double probabilityParameter = 1.0/(problem->getNumberOfVariables());
  double distributionIndexParameter = 20.0;
  parameters["probability"] =  &probabilityParameter;
  parameters["distributionIndex"] = &distributionIndexParameter;
  mutation = new PolynomialMutation(parameters);
  // Add the operators to the algorithm
  algorithm->addOperator("mutation",mutation);

  // Add the indicator object to the algorithm
  //algorithm->setInputParameter("indicators", indicators) ;

  // Execute the Algorithm
  t_ini = clock();
  SolutionSet * population = algorithm->execute(); 
  t_fin = clock();
  double secs = (double) (t_fin - t_ini);
  secs = secs / CLOCKS_PER_SEC;

  // Result messages
  cout << "Total execution time: " << secs << "s" << endl;
  cout << "Variables values have been written to file VAR" << endl;
  population->printVariablesToFile("VAR");
  cout << "Objectives values have been written to file FUN" << endl;
  population->printObjectivesToFile("FUN");

  delete mutation;
  delete population;
  delete algorithm;

} // main


