#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cfloat>
using namespace std;

const double EULER_CONSTANT = exp(1.0);

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

const int k = 5;

double phi1 = 2.05;
double phi2 = 2.05;

string top;
string problem;
int nSize;

string medians = "";
std::vector<std::vector<double> > bestValuesForTestVectorOfVectors;

struct particle{
    vector<double> position;
    vector<double> velocity;
    
    double bestValue = DBL_MAX;
    vector<double> bestPosition;
    
    double bestNeighborhoodValue = DBL_MAX;
    vector<double> bestNeighborhoodPosition;
};

vector<particle> partVect;

void genParticleVector();
void genNeighborhoods();
void performIteration();
void getNBest();
void genVN();

void martinTests();

double evaluateAckley(particle particleStruct);
double evaluateRok(particle particleStruct);
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

vector<double> evaluateSphereVector(vector<particle> particleVector);
double evaluateSphereVector(particle z);
double getMedian(vector<double> testResults);

int main(int argc, char* argv[]){
    
    srand(time(0));
//
//    bestValueFound = DBL_MAX;
//    
//    top = argv[1];
//    numParticles = atoi(argv[2]);
//    numIterations = atoi(argv[3]);
//    problem = argv[4];
//    numDimensions = atoi(argv[5]);
//    
//    setBounds();
//    setTop();
//    
//    genParticleVector();
//    
//    if(top != "gl"){		//if the topology type is global, no need for neighborhoods
//        genNeighborhoods();
//    }
//    
//    for(int i = 0; i < numParticles; i++){	//initializes "neighborhoodBest" (nBest) vector
//        nBest.push_back(-1);
//    }
//    
//    runPSO();
//    cout<<bestValueFound<<endl;

      martinTests();
    
}

void martinTests(){
    
    vector <string> problems;
    vector <string> topologies;
    vector <int> particleCounts;
    
    problems.push_back("ras");
   // problems.push_back("rok");
    //problems.push_back("sph");
   // problems.push_back("ack");


    topologies.push_back("ri");
    topologies.push_back("vn");
    topologies.push_back("ra");
    topologies.push_back("gl");
    
    numDimensions = 30;
    numIterations = 10000;


    
   particleCounts.push_back(16);
   particleCounts.push_back(30);
   particleCounts.push_back(49);
  //particleCounts.push_back(50);
    
    // how many times we want to run each one
    int desiredTests = 20;
    
    for (int i = 0; i < problems.size(); i++) {
        
        for (int j = 0; j < topologies.size(); j++) {
            
            for (int m = 0; m < particleCounts.size(); m++) {
                
                vector<double> testResults;
                
                for (int n = 0; n < desiredTests; n++) {
                  //  cout<<n<<endl;
                    bestValueFound = DBL_MAX;
                    
                    problem = problems.at(i);
                    numParticles = particleCounts.at(m);
                    top = topologies.at(j);
                    
                    setBounds();	//sets upper and lower position and velocity bounds specific to each problem
                    setTop();		//sets the neighborhood size var depending on which neighborhood type is being used
                    
                    
                    genParticleVector();	//generates the particle vector
                    
                    if(top != "gl"){		//if the topology type is global, no need for neighborhoods
                        genNeighborhoods();
                    }
                    
                    
                    for(int i = 0; i < numParticles; i++){	//initializes "neighborhoodBest" (nBest) vector
                        nBest.push_back(-1);
                    }
                    
                    runPSO();

                    testResults.push_back(bestValueFound);
                    
                    nBest.clear();
                    partVect.clear();
                    //empty particle vector
                    //empty nBest
                }
                
                double medianResult;
                medianResult = getMedian(testResults);
                // for(int i = 0; i < testResults.size(); i++){
                //     cout<<testResults[i]<<endl;
                // }
                
                double averageResult = 0;
                
                for(int k = 0; k < testResults.size(); k++){
                    averageResult += testResults[k];
                }
                averageResult = averageResult / testResults.size();

                std::vector<double> mediansForTenTests;



                std::vector<double> mediansForPrinting;
                for(int k = 0; k < bestValuesForTestVectorOfVectors[0].size(); k++){
                    for(int z = 0; z < bestValuesForTestVectorOfVectors.size(); z++){
                        mediansForTenTests.push_back(bestValuesForTestVectorOfVectors[z][k]);
                    }
                    mediansForPrinting.push_back(getMedian(mediansForTenTests));
                    mediansForTenTests.clear();
                }

                // You can easily change the format of the print outs here if you
                // want to change it up and make excel easier
                cout<<problem<<","<<top<<","<<numParticles<<","<<numDimensions<<","<<numIterations<<","<<averageResult<<","<<medianResult<<","<<
                mediansForPrinting[0]<<","<<mediansForPrinting[1]<<","<<mediansForPrinting[2]<<","<<mediansForPrinting[3]<<","
                <<mediansForPrinting[4]<<","<<mediansForPrinting[5]<<","<<mediansForPrinting[6]<<","<<mediansForPrinting[7]<<","
                <<mediansForPrinting[8]<<","<<mediansForPrinting[9]<<","<<mediansForPrinting[10]<<endl;

                // cout << "problem: " << problem;
                // cout << ", topology: " << top;
                // cout << ", numParticles: " << numParticles;
                // cout << " --- median best value of " << desiredTests << " tests = " << medianResult << endl;
                // cout << " --- bestValueFound: "<<bestValueFound<<endl;
                // cout << " --- average best value of " <<desiredTests << " tests = " << averageResult << endl;

                bestValuesForTestVectorOfVectors.clear();


            }
        }
    }
    
}

double getMedian(vector<double> testResults) {
    
    sort(testResults.begin(), testResults.end());
    int medIndex = testResults.size() / 2;
    return testResults.at(medIndex);
}


void runPSO(){
    std::vector<double> bestValuesForTest;

    for(int i = 0; i < numIterations; i++){
        performIteration();
        if(top == "ra"){		//extra step when using random topology
            reRandomizeNeighborhoods();
        }
        if(i % 1000 == 0){
            bestValuesForTest.push_back(bestValueFound);
        }
    }
    bestValuesForTest.push_back(bestValueFound);

    bestValuesForTestVectorOfVectors.push_back(bestValuesForTest);
}

void genParticleVector(){
    double posVal;
    double velVal;
    
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

void genNeighborhoods(){		//wrapper function
    if(top == "ri"){			//ring is short so I just wrote it here
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

void performIteration(){
    //prep for update
    evaluations = getErrors(partVect, problem);		//evaluations is a vector that stores to calculated value of each particle's position
    gBest = getIndexOfMinParticle(evaluations);		//gBest is the index of the global best particle
    
    getNBest();                                     //sets the nBest vector for all cases
    
    
    for(int i = 0; i < numParticles; i++){
        if(evaluations[i] < partVect[i].bestValue){
            partVect[i].bestValue = evaluations[i];
            partVect[i].bestPosition = partVect[i].position;
        }
    }
    
    
    
    
    if(evaluations[gBest] < bestValueFound){
        bestValueFound = evaluations[gBest];			//keeps track of best value discovered so far
    }
    
    //cout << endl << bestValueFound;
    
    vector<particle> partVectCopy = partVect;			//creates a copy of the particle vector so that all particles can be updated before the update "takes effect"
    //actual update
    for(int i = 0; i < numParticles; i++){
        for(int j = 0; j < numDimensions; j++){
            
            //update formulas in paper from Majercik
            
            // double firstPart = (phi1 * (double)rand()/(double)RAND_MAX) * (partVectCopy[nBest[i]].position[j] - partVectCopy[i].position[j]);
            // double secondPart = (phi2 * (double)rand() / (double)RAND_MAX * (partVectCopy[gBest].position[j] - partVectCopy[i].position[j]));
            // double result = 0.7298 * (partVectCopy[i].velocity[j] + firstPart + secondPart);
            // partVect[i].velocity[j] = result;
            
            partVect[i].velocity[j] = 0.7298 * (partVectCopy[i].velocity[j] + (phi1 * ((double)rand() / (double)RAND_MAX) * (partVectCopy[i].bestPosition[j] - partVectCopy[i].position[j]))
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


void getNBest(){			//sets the neighborhood best (nBest) vector
    vector<int> temp;
    if(top  == "gl"){
        for(int i = 0; i < numParticles; i++){
            temp.push_back(gBest);
            if(evaluations[gBest] < partVect[i].bestNeighborhoodValue){
                partVect[i].bestNeighborhoodValue = evaluations[gBest];
                partVect[i].bestNeighborhoodPosition = partVect[gBest].position;
            }
        }
        nBest = temp;
    }else{
        //else case
        int neighborhoodBest;
        for(int i = 0; i < numParticles; i++){
            neighborhoodBest = neighborhoods.at(i).at(0);
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

void setTop(){
    if(top == "vn" || top == "ra"){
        nSize = 5;
    }else if(top == "ri"){
        nSize = 3;
    }
}
//JAMES CODE BELOW

std::vector<double> getErrors(std::vector<particle> particleVector, std::string problemName){
    
    if(problemName == "rok"){
        return evaluateRokVector(particleVector);
    }
    else if(problemName == "ack"){
        return evaluateAckleyVector(particleVector);
    }else if(problemName == "sph"){
        return evaluateSphereVector(particleVector);
    }
    return evaluateRastriginVector(particleVector);
}

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

//take a single particle, spit out double
//for evaluating for each function

//take a vector of particles and return index of best particle

//position is position in each dimensions
//veolicy is veolictty in each dimensions

//works with two particles with x = 8, y = 8


std::vector<double> evaluateRokVector(std::vector<particle> particleVector){
    std::vector<double> vectorOfEvaluations;
    
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


//evaluates correctly for two particles with x=8, y=8
std::vector<double> evaluateAckleyVector(std::vector<particle> particleVector){
    
    double PI = M_PI;
    std::vector<double> vectorOfEvaluations;
    
    for(int z = 0; z < particleVector.size(); z++){
        vectorOfEvaluations.push_back(evaluateAckley(particleVector[z]));
        
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

//STATUS: should evaluate to 4 at x=2, and it does
std::vector<double> evaluateRastriginVector(std::vector<particle> particleVector){
    
    
    std::vector<double> vectorOfEvaluations;
    
    for(int z = 0; z < particleVector.size(); z++){
        vectorOfEvaluations.push_back(evaluateRastrigin(particleVector[z]));
    }
    
    return vectorOfEvaluations;
}

double evaluateSphere(particle z){
    double sum = 0;
    for(int i = 0; i < numDimensions; i++){
        sum = sum + z.position[i] * z.position[i];
    }
    return sum;
}

vector<double> evaluateSphereVector(vector<particle> particleVector){
    vector<double> z;
    for(int i = 0; i < numParticles; i++){
        z.push_back(evaluateSphere(particleVector[i]));
    }
    return z;
}


//MARTIN CODE BELOW

//these functions should change the neighborhoods vector of vector of ints so that each index of neighborhoods
//represents that particle's neighborhood. (i.e., if the vector of ints at index 0 of neighborhoods is equal to {0, 3, 56, 2}, then
//we know that the particles indexed at 3, 56, and 2 in "partVect" are neighbors of the particle indexed at 0 in "partVect").
//Each neighborhood should contain its own index (i.e., the vector of ints at index 0 of neighborhoods should contain 0.)


// populates the global var neighborhoods with random neighborhoods (no dupes in a neighborhood)
void genRandomNeighborhoods(){
    
    neighborhoods.clear();
    for(int i = 0; i < numParticles; i++){
        vector<int> neighborhood;
        neighborhood.push_back(i);
        
        while (neighborhood.size() < k) {
            int randomParticle = rand() % numParticles;
            if (find(neighborhood.begin(), neighborhood.end(), randomParticle) == neighborhood.end()) {
                neighborhood.push_back(randomParticle);
            }
        }
        neighborhoods.push_back(neighborhood);
    }
}

// call this on each iteration. It rerandomizes 20% of the random neighborhoods that already exist (from handout)
// theres a chance i interpreted the handout wrong on this one because the wording is a little confusing
// but i think its right
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


// populates Neighbors with a VN topology. The special cases get a bit messy, but work.
// some of the edge cases end up with duplicate neighbors, but this is fine and unavoidable really.
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

double evaluateRok(particle particleStruct){
    
    double evaluation = 0;
    
    for(int i = 0; i < particleStruct.position.size()-1; i++){
        double positionAtI = particleStruct.position[i];
        double product = particleStruct.position[i+1] - positionAtI * positionAtI;
        
        evaluation += (100 * (product * product) + (positionAtI -1) * (positionAtI - 1));
    }
    
    return evaluation;
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
    //std::cout<<firstBracket<<std::endl;
    
    for(int i = 0; i < particleStruct.position.size(); i++){
        secondBracket += cos(2 * PI * particleStruct.position[i]);
    }
    secondBracket *= (1.00f / (double) particleStruct.position.size());
    //std::cout<<secondBracket<<std::endl;
    
    finalSum = -20.00f * pow(EULER_CONSTANT, firstBracket) - pow(EULER_CONSTANT, secondBracket) + 20.00f + EULER_CONSTANT;
    
    return finalSum;
}
