#include "Data.hpp"
#include <iostream>
#include <vector>

using namespace std;

typedef struct {
    vector<int> sequence;
    double cost;
} TspSolution;

struct InsertionInfo {
    //Node that'll be inserted as graph gets built
    int insertedNode;
    //Edge that'll be removed as insertedNode gets added into the graph
    //To be substituted by 2 other edges. Both of which connect one previously-interconnected
    //node to to the new insertedNode 
    int removedGraphEdge;
    //Cost of inserting this new node into the graph, among tho other nodes
    double cost;
};

//BuildSolution Utility funcs
/**
 * @brief Calcs the cost of inserting the new node into the graph. (vertex and node used interchangaebly here)
 * @param tspSol To-be-evaluated TSP Solution
 * @param unaddedVertices Vertices / Nodes that still don't belong to the solution subtour
 * @param distMatrix Matrix holding distances among graph nodes
 * @return insertionCost Cost of inserting nodes into TSP Solution graph
 */
vector<InsertionInfo> CalcNodeInsertionCost(TspSolution& tspSol, vector<int>& unaddedVertices, double** distMatrix)
{
    //Specifies size of vec
    vector<InsertionInfo> insertionCost = vector<InsertionInfo>((tspSol.sequence.size() - 1) * unaddedVertices.size());
    int l = 0;
    int i, j;

    //Iterating among all nodes in the solution sequence, going edge-by-edge (df 2 nodes on each edge)
    for(int edgeIter = 0; edgeIter < (int)tspSol.sequence.size() - 1; edgeIter++){
        i = tspSol.sequence[edgeIter];  //Vertex 1
        j = tspSol.sequence[edgeIter + 1];  //Vertex 2

        //For each node/vertex
        for(int k : unaddedVertices){
            insertionCost[l].cost = distMatrix[i][k] + distMatrix[j][k] - distMatrix[i][j];    //TODO: Check about c
            insertionCost[l].insertedNode = k;
            insertionCost[l].removedGraphEdge = edgeIter;
            l++;
        }
    }

    return insertionCost;
}

//Chooses 3 random nodes in the graph from 1 up to dimension (numOfNodes)
//Makes a sequence from them following {1, x, y, z, 1} and puts the 3 added nodes into solAddedNodes
//solAddedNodes must be empty
vector<int> Choose3RandNodes(int dimension, vector<int>& solAddedNodes)
{
    vector<int> starterSeq;
    int randNode = 0;
    bool nodeAlreadyAddedFlag = false;

    srand(static_cast<unsigned int>(time(0)));

    starterSeq.reserve(5);  //{1, x, y, z, 1}
    starterSeq.push_back(1);

    //Choosing ramdonly 'till 3 non-equal nodes are added
    //Doing this cause i can't add the same node twice (leads to INFINITE distance)
    while((int)starterSeq.size() != 4){
        //Randomly chooses a node from 0 up to numOfNodes
        randNode = rand() % (dimension + 1);

        //Handles if rand chooses 0 or 1, since distances 1-indexed here and 1 is the first and last cities
        //TODO: Check if it can be done on the randNode attr line straight up
        if((randNode == 0) || (randNode == 1)){
            randNode = 2;
        }

        //Checks if this node has already been added to starterSeq
        for(int i = 0; i < (int)solAddedNodes.size(); i++){
            if(randNode == solAddedNodes[i]){
                nodeAlreadyAddedFlag = true;
                break;
            }
        }

        //Only adds it if this node isn't a repeat
        if(!nodeAlreadyAddedFlag){
            starterSeq.push_back(randNode);
            solAddedNodes.push_back(randNode);
        }
        nodeAlreadyAddedFlag = false; //For next iter
    }
    starterSeq.push_back(1);

    return starterSeq;
}

//Receives a vector with the ramdomly-added-to-solution nodes, and removes them from the vector containing all nodes in the graph
//Effectively returns CL
vector<int> GetUnchosenNodes(int dimension, vector<int>& solAddedNodes)
{
    vector<int> unaddedNodes;
    bool nodeFound;

    //Iterates through s from 2 to numOfNodes (distMatrix is 1-indexed and 1 will always be in the solution)
    for(int i = 2; i <= dimension; i++){
        //Checks if current array element is any of the 3 in solAddedNodes (except for 1, since it's beginning/end)
        for(int j = 0; j < (int)solAddedNodes.size(); j++){
            //Check if this node was added to solSequence, immediately breaking out of the 3-check loop if it's found
            if(i == solAddedNodes[j]){
                nodeFound = true;
                break;
            }
            nodeFound = false;
        }

        if(!nodeFound){
            unaddedNodes.push_back(i);
        }
    }
    return unaddedNodes;
}

/**
 * @brief Builds a fair solution (though still far from optimized) Using Greedy Randomized Adaptive Search (Insertion-by-cheapest)
 * @return TspSolution
 */
TspSolution BuildSolution(double** distMatrix, int dimension)
{
    TspSolution tspSol;
    vector<int> addedNodes;
    tspSol.sequence = Choose3RandNodes(dimension, addedNodes);   //gets s
    vector<int> unaddedNodes = GetUnchosenNodes(dimension, addedNodes);   //gets CL

    while(!unaddedNodes.empty()){
        vector<InsertionInfo> insertionCost = CalcNodeInsertionCost(tspSol, unaddedNodes, distMatrix);

        //insertionCost;

        double alpha = (double)rand() / RAND_MAX;

        int selected = rand() % ((int)ceil(alpha * insertionCost.size()));

        //insertIntoSolution(tspSol, insertionCost[selected].insertedNode); //k
    }

    return tspSol;
}

TspSolution IteratedLocalSearch(int maxIters, int maxIterILS, Data& data)
{
    TspSolution bestOfAllSolution;
    bestOfAllSolution.cost = INFINITY;

    for(int i = 0; i < maxIters; i++){
        //Builds a beginning solution based on fair guesses
        TspSolution currIterSolution = BuildSolution(data.getMatrixCost(), data.getDimension());
        TspSolution currBestSolution = currIterSolution;

        int iterILS = 0;

        while(iterILS <= maxIterILS){
            //Tries to enhance the fair-guessed solution
            //By doing small modifications to it
            //LocalSearch(&currIterSolution);

            if(currIterSolution.cost < currBestSolution.cost){
                currBestSolution = currIterSolution;
                iterILS = 0;
            }

            //If not possible to make it better, shake the current best solution up a lil
            //to see if we didn't just go into a 'local best pitfall'
            //currIterSolution = Disturbance(currBestSolution);
            iterILS++;
        }

        if(currBestSolution.cost < bestOfAllSolution.cost){
            bestOfAllSolution = currBestSolution;
        }
    }

    return bestOfAllSolution;
}

int main(int argc, char** argv)
{
    auto data = Data(argc, argv[1]);
    data.read();

    size_t n = data.getDimension();

    cout << "Dimension: " << n << endl;
    cout << "DistanceMatrix: " << endl;
    data.printMatrixDist();

    //Heuristic goes here
    IteratedLocalSearch(5, 5, data);

    cout << "----------------------\n";
    cout << "Exemplo de Solucao s = ";

    double cost = 0.0;

    for(size_t i = 1; i < n; i++){
        cout << i << " -> ";
        cost += data.getDistance(i, i+1);
    }

    cost += data.getDistance(n, 1);

    cout << n << " -> " << 1 << endl;
    cout << "Custo de S: " << cost << endl;

    return 0;
}
