#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cfloat>
#include <time.h>
using namespace std;

const double EULER_CONSTANT = exp(1.0);


// Global Variables
int numDimensions;
int numParticles;
int numIterations;
vector< vector<int> > neighborhoods;
vector<double> evaluations;
vector<int> nBest;
int gBest;
double bestValueFound = DBL_MAX;

double lpbound;
double upbound;
double lvbound;
double uvbound;

const int k = 5; // size of random neighborhood

double phi1 = 2.05;
double phi2 = 2.05;

string top; //topology as string
string problem; //problem as string
int nSize;  //neighborhood size

// particle stuct houses current position and velocity, as well as neighborhood and personal bests
struct particle{
    vector<double> position;
    vector<double> velocity;
    
    double bestValue = DBL_MAX;
    vector<double> bestPosition;
    
    double bestNeighborhoodValue = DBL_MAX;
    vector<double> bestNeighborhoodPosition;
};

vector<particle> partVect;

// FUNCTION DECLARATIONS
void genParticleVector();
void genNeighborhoods();
void performIteration();
void getNBest();

double evaluateAckley(particle particleStruct);
std::vector<double> evaluateRastriginVector(std::vector<particle> particleVector);
std::vector<double> evaluateRokVector(std::vector<particle> particleVector);
std::vector<double> evaluateAckleyVector(std::vector<particle> particleVector);
double evaluateRastrigin(particle newParticle);


std::vector<double> getErrors(std::vector<particle> particleVector, std::string problemName);
int getIndexOfMinParticle(std::vector<double> vectorOfEvaluations);

void genVN();
void genRandomNeighborhoods();
void reRandomizeNeighborhoods();
void printNeighborhoods();
void setBounds();
void setTop();

void runPSO();



// -------------------------------------------------------------------------------


int main(int argc, char* argv[]){

    if(argc < 6){
        cout<< "Usage: ./a.out topology numberOfParticles numberOfIterations problemType numberOfDimensions" << endl;
        exit(-1);
    }
    
    srand(time(0));
    
    bestValueFound = DBL_MAX;

    //getting arguments from command line
    
    top = argv[1];
    numParticles = atoi(argv[2]);
    numIterations = atoi(argv[3]);
    problem = argv[4];
    numDimensions = atoi(argv[5]);

    clock_t start = clock();    

    //prints for problem
    cout << "Topology: " << top << endl;
    cout << "Swarm size: " << numParticles << endl;
    cout << "Number of iterations: " << numIterations << endl;
    cout << "Problem function: " << problem <<endl;
    cout << "Number of dimensions: " << numDimensions <<endl;

    // set bounds and topology based on command line
    setBounds();
    setTop();
    genParticleVector();
    if(top != "gl"){		//if the topology type is global, no need for neighborhoods
        genNeighborhoods();
    }

    for(int i = 0; i < numParticles; i++){	//initializes "neighborhoodBest" (nBest) vector
        nBest.push_back(-1);
    }

    runPSO();

    clock_t stop = clock();    
    double elapsed = (double)(stop - start) * 1000.0f / CLOCKS_PER_SEC / 1000.0f;  

    cout << "Elapsed time: " << elapsed << " seconds" << endl;
    cout << endl;
    cout << "Best value found: " << bestValueFound << endl;
    cout << endl;
    for(int i = 0; i < evaluations.size(); i++){
        cout<<"Final value found for particle " << i << " = " <<evaluations[i] << endl;
    }    
}


//generates a vector of particles with randomized position and velocity within the bounds
void genParticleVector(){
    
    double posVal;
    double velVal;
    
    // uses specific bounds from the problem to initialize the particles
    for(int i = 0; i < numParticles; i++){
        particle newParticle;
        for(int j = 0; j < numDimensions; j++){
            posVal = (upbound - lpbound) * (double)rand()/(double)RAND_MAX + lpbound;	//creates a random position value within the position bounds and adds it to the position vector
            newParticle.position.push_back(posVal);
            velVal = (uvbound - lvbound) * (double)rand()/(double)RAND_MAX + lvbound;	//creates a random velocity value within the velocity bounds and adds it to the velocity vector
            newParticle.velocity.push_back(velVal);
        }
        newParticle.bestPosition = newParticle.position;
        partVect.push_back(newParticle);
    }
}

void genNeighborhoods(){		//neighborhoodwrapper function
    
    // horizontal neighbors and self
    if(top == "ri"){			
        for(int i = 0; i < numParticles; i++){
            vector<int> neighborhood;
            if(i==0){
                neighborhood.push_back(i);
                neighborhood.push_back(i + 1);
                neighborhood.push_back(numParticles - 1);
            }else if(i == numParticles - 1){
                neighborhood.push_back(i - 1);
                neighborhood.push_back(i);
                neighborhood.push_back(0);
            }else{
                neighborhood.push_back(i - 1);
                neighborhood.push_back(i);
                neighborhood.push_back(i + 1);
            }
            neighborhoods.push_back(neighborhood);
        }
    }else if(top == "vn"){
        genVN();
    }else if(top == "ra"){
        genRandomNeighborhoods();
    }else if(top == "gl"){
        //special case
    }else{cout << "ERROR";}
}


// -------------------------------------------------------------------------------

// -------------------------------------------------------------------------------


void runPSO(){
    for(int i = 0; i < numIterations; i++){
        performIteration();
        if(top == "ra"){        //extra step when using random topology
            reRandomizeNeighborhoods();
        }
    }
}

//wrapper function for a single iteration of PSO
//includes getting function values, getting the best particle, gets neighborhood best, and then updates
void performIteration(){
    //prep for update
    evaluations = getErrors(partVect, problem);		//evaluations is a vector that stores to calculated value of each particle's position
    gBest = getIndexOfMinParticle(evaluations);		//gBest is the index of the global best particle
    getNBest();                                     //sets the nBest vector for all cases
    double CONSTRICTION_FACTOR = 0.7298; //from handout, never changes because phi values are constant

    //update personal best value and position for each particle
    for(int i = 0; i < numParticles; i++){
        if(evaluations[i] < partVect[i].bestValue){
            partVect[i].bestValue = evaluations[i];
            partVect[i].bestPosition = partVect[i].position;
        }
    }

    if(evaluations[gBest] < bestValueFound){
        bestValueFound = evaluations[gBest];	//keeps track of best value discovered so far
    }    
    vector<particle> partVectCopy = partVect;	//creates a copy of the particle vector so that all particles can be updated before the update "takes effect"
    //actual update of position
    for(int i = 0; i < numParticles; i++){
        for(int j = 0; j < numDimensions; j++){
            
            //update formula
            partVect[i].velocity[j] = CONSTRICTION_FACTOR * (partVectCopy[i].velocity[j] + (phi1 * ((double)rand() / (double)RAND_MAX) * (partVectCopy[i].bestPosition[j] - partVectCopy[i].position[j]))
                                                + (phi2 * ((double)rand() / (double)RAND_MAX) * (partVectCopy[i].bestNeighborhoodPosition[j] - partVectCopy[i].position[j])));
            //keeping the velocity bounded
            if(partVect[i].velocity[j] > uvbound){
                partVect[i].velocity[j] = uvbound;
            }else if(partVect[i].velocity[j] < lvbound){
                partVect[i].velocity[j] = lvbound;
            }

            partVect[i].position[j] = partVect[i].position[j] + partVect[i].velocity[j];
        }
    }
}


//finds the neighborhood best for global case and all other cases.
void getNBest(){

    vector<int> temp;
    
    if(top  == "gl"){ // neighborhood best is just the global best value we already got
        for(int i = 0; i < numParticles; i++){
            temp.push_back(gBest);
            if(evaluations[gBest] < partVect[i].bestNeighborhoodValue){
                partVect[i].bestNeighborhoodValue = evaluations[gBest];
                partVect[i].bestNeighborhoodPosition = partVect[gBest].position;
            }
        }
        nBest = temp;
    }else{
        //else case, goes through each neighborhood and finds its best position
        int neighborhoodBest;
        for(int i = 0; i < numParticles; i++){
            neighborhoodBest = neighborhoods.at(i).at(0);
            // update each particle and updates its neighborhood best value
            for(int j = 0; j < nSize; j++){
                double currMin = evaluations[neighborhoodBest];
                if(evaluations[neighborhoods.at(i).at(j)] < currMin){
                    neighborhoodBest = neighborhoods.at(i).at(j);
                    currMin = evaluations[neighborhoodBest];
                }
            }
            temp.push_back(neighborhoodBest);
            if(evaluations[neighborhoodBest] < partVect[i].bestNeighborhoodValue){
                partVect[i].bestNeighborhoodValue = evaluations[neighborhoodBest];
                partVect[i].bestNeighborhoodPosition = partVect[neighborhoodBest].position;
            }
        }
        nBest = temp;
    }
}


//returns index of particle with best value
int getIndexOfMinParticle(std::vector<double> vectorOfEvaluations){
    
    int bestIndex = 0;
    //arbitrarily set
    double minParticleValue = vectorOfEvaluations[bestIndex];
    
    for(int i = 0; i < vectorOfEvaluations.size(); i++){
        if(vectorOfEvaluations[i] < minParticleValue){
            minParticleValue = vectorOfEvaluations[i];
            bestIndex = i;
        }
    }
    return bestIndex;
}

// -------------------------------------------------------------------------------
// FUNTION EVALUATIONS
// -------------------------------------------------------------------------------

//function is a wrapper for evaluating the right problem function
std::vector<double> getErrors(std::vector<particle> particleVector, std::string problemName){
    
    if(problemName == "rok"){
        return evaluateRokVector(particleVector);
    }
    else if(problemName == "ack"){
        return evaluateAckleyVector(particleVector);
    }

    return evaluateRastriginVector(particleVector);
}

//returns a vector of evaluations for rosenbrock
std::vector<double> evaluateRokVector(std::vector<particle> particleVector){
    std::vector<double> vectorOfEvaluations;
    
    // for each particle
    for(int z = 0; z < particleVector.size(); z++){
        double evaluation = 0;
        
        for(int i = 0; i < particleVector[z].position.size()-1; i++){

            double product = particleVector[z].position[i+1] - particleVector[z].position[i] * particleVector[z].position[i];
            
            evaluation += (100 * (product * product) + (particleVector[z].position[i] -1) * (particleVector[z].position[i] - 1));
        }
        vectorOfEvaluations.push_back(evaluation);
    }
    
    return vectorOfEvaluations;
}

std::vector<double> evaluateAckleyVector(std::vector<particle> particleVector){
    
    double PI = M_PI;
    std::vector<double> vectorOfEvaluations;
    
    for(int z = 0; z < particleVector.size(); z++){
        vectorOfEvaluations.push_back(evaluateAckley(particleVector[z]));
        
    }
    
    return vectorOfEvaluations;
}

std::vector<double> evaluateRastriginVector(std::vector<particle> particleVector){
    
    
    std::vector<double> vectorOfEvaluations;
    
    for(int z = 0; z < particleVector.size(); z++){
        vectorOfEvaluations.push_back(evaluateRastrigin(particleVector[z]));
    }
    
    return vectorOfEvaluations;
}

double evaluateRastrigin(particle newParticle){
    
    double PI = M_PI;

    double sum = 0;
    
    for(int i = 0; i < newParticle.position.size(); i++){
        double positionAtI = newParticle.position[i];
        sum+= ( (positionAtI * positionAtI) - (10 * cos(2*PI*positionAtI)) + 10) ;
    }
    return sum;
}

double evaluateAckley(particle particleStruct){
    
    double evaluation = 0;
    
    double firstBracket = 0;
    double secondBracket=0;
    
    double finalSum = 0;
    
    
    double PI = M_PI;
    
    for(int i = 0; i < particleStruct.position.size(); i++){
        firstBracket += (particleStruct.position[i] * particleStruct.position[i]);
    }
    firstBracket *= (1.00f / (double) particleStruct.position.size() );
    firstBracket = sqrt(firstBracket);
    firstBracket *= -0.2;
    
    for(int i = 0; i < particleStruct.position.size(); i++){
        secondBracket += cos(2 * PI * particleStruct.position[i]);
    }
    secondBracket *= (1.00f / (double) particleStruct.position.size());
    
    finalSum = -20.00f * pow(EULER_CONSTANT, firstBracket) - pow(EULER_CONSTANT, secondBracket) + 20.00f + EULER_CONSTANT;
    
    return finalSum;
}

// -------------------------------------------------------------------------------
// END EVALUATIONS
// -------------------------------------------------------------------------------

// -------------------------------------------------------------------------------
// GEN NEIGHBORHOODS
// -------------------------------------------------------------------------------

// populates the global var neighborhoods with random neighborhoods (no dupes in a neighborhood)
void genRandomNeighborhoods(){
    
    neighborhoods.clear();
    for(int i = 0; i < numParticles; i++){
        vector<int> neighborhood;
        neighborhood.push_back(i);
        
        // k = size of random neighborhood
        while (neighborhood.size() < k) {
            int randomParticle = rand() % numParticles;
            // check for dupes
            if (find(neighborhood.begin(), neighborhood.end(), randomParticle) == neighborhood.end()) {
                neighborhood.push_back(randomParticle);
            }
        }
        neighborhoods.push_back(neighborhood);
    }
}

// call this on each iteration. It rerandomizes 20% of the random neighborhoods that already exist (from handout)

void reRandomizeNeighborhoods() {
    
    for(int i = 0; i < numParticles; i++){
        int randomProb = rand() % 100;
        if(randomProb < 20) {
            neighborhoods[i].clear();
            neighborhoods[i].push_back(i);
            
            while (neighborhoods[i].size() < k) {
                int randomParticle = rand() % numParticles;
                if (find(neighborhoods[i].begin(), neighborhoods[i].end(), randomParticle) == neighborhoods[i].end()) {
                    neighborhoods[i].push_back(randomParticle);
                }
            }
        }
    }
}


// populates Neighbors with a VN topology. The special cases get a bit messy, but it works for populations of any size
void genVN(){
    
    neighborhoods.clear();
    
    //gets dimensions of 2d array
    int root = pow(numParticles, 0.5);
    if (root * root != numParticles) {
        root++;
    }
    
    //makes array with numParticles entries (all other entries are -1)
    int array[root][root];
    int particlesLeft = numParticles;

    for(int col = 0; col < root; col++){
        for(int row = 0; row < root; row++){
            if (particlesLeft > 0) {
                array[col][row] = col * root + row;
                particlesLeft--;
            } else {
                array[col][row] = -1;
            }
        }
    }
    
    //puts empty neighborhoods in Neighborhoods
    for(int i = 0; i < numParticles; i++){
        vector<int> neighborhood;
        neighborhoods.push_back(neighborhood);
    }
    
    for(int col = 0; col < root; col++){
        for(int row = 0; row < root; row++){
            //make sure we are on a particle
            if (array[col][row] != -1) {
                //add particle to its own neighborhood
                neighborhoods[array[col][row]].push_back(array[col][row]);
                
                //up of particle
                if(col == 0){
                    if (array[root - 1][row] != -1) {
                        neighborhoods[array[col][row]].push_back(array[root - 1][row]);
                    } else {
                        int offset = 2;
                        while (array[root - offset][row] == -1) {
                            offset++;
                        }
                        neighborhoods[array[col][row]].push_back(array[root - offset][row]);
                    }
                }else{
                    neighborhoods[array[col][row]].push_back(array[col - 1][row]);
                }
                
                //down of particle
                if(col == root - 1){
                    neighborhoods[array[col][row]].push_back(array[0][row]);
                }else{
                    if (array[col + 1][row] != -1) {
                        neighborhoods[array[col][row]].push_back(array[col + 1][row]);
                    } else {
                        neighborhoods[array[col][row]].push_back(array[0][row]);
                    }
                }
                
                //left of particle
                if(row == 0){
                    if (array[col][root - 1] != -1) {
                        neighborhoods[array[col][row]].push_back(array[col][root - 1]);
                    } else {
                        int offset = 2;
                        while (array[col][root - offset] == -1) {
                            offset++;
                        }
                        neighborhoods[array[col][row]].push_back(array[col][root - offset]);
                    }
                }else{
                    neighborhoods[array[col][row]].push_back(array[col][row - 1]);
                }
                
                //right of particle
                if(row == root - 1){
                    neighborhoods[array[col][row]].push_back(array[col][0]);
                }else{
                    if (array[col][row + 1] != -1) {
                        neighborhoods[array[col][row]].push_back(array[col][row + 1]);
                    } else {
                        neighborhoods[array[col][row]].push_back(array[col][0]);
                    }
                }
            }
        }
    }   
}

// function to print neighborhoods list
void printNeighborhoods() {
    
    cout << "list of neighborhoods:" << endl;
    for (int i = 0; i < numParticles; i++) {
        cout << "neighborhood for particle " << i << ": ";
        for (int j = 0; j < nSize; j++) {
            cout << neighborhoods.at(i).at(j) << ",";
        }
        cout << endl;
    }
}


//sets bounds depending upon function
void setBounds(){
    if(problem == "rok"){
        upbound = 30.0;
        lpbound = 15.0;
        uvbound = 2.0;
        lvbound = -2.0;
    }else if(problem == "ack"){
        upbound = 32.0;
        lpbound = 16.0;
        uvbound = 4.0;
        lvbound = -2.0;
    }else if(problem == "ras"){
        upbound = 5.12;
        lpbound = 2.56;
        uvbound = 4.0;
        lvbound = -2.0;
    }else{
        upbound = 15;
        lpbound = 15;
        uvbound = 2.0;
        lvbound = -2.0;
    }
}

//sets neighborhood size for topology
void setTop(){
    if(top == "vn" || top == "ra"){
        nSize = 5;
    }else if(top == "ri"){
        nSize = 3;
    }
}