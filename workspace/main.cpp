#include "Data.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

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
        //Checks if current array element is any of the 3 in solAddedNodes
        for(int j = 0; j < (int)solAddedNodes.size(); j++){
            //Checks if this node was added to solSequence, immediately breaking out of the 3-check loop if it's found
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

//Sorting the costs
bool CompareByCost(InsertionInfo& prevInsertionInfo, InsertionInfo& currInsertionInfo)
{
    return prevInsertionInfo.cost < currInsertionInfo.cost;
}

void SortAscendingByCost(vector<InsertionInfo>& insertionInfo)
{
    sort(insertionInfo.begin(), insertionInfo.end(), CompareByCost);
}

//Inserts a node into the tsp solution
//And sums the current tspSol cost with the inserted node's own
void InsertIntoSolution(TspSolution& tspSol, InsertionInfo& nodeInsertionInfo)
{
    //First adds up the cost of the to-be-added node
    tspSol.cost += nodeInsertionInfo.cost;

    //Adds the node into the tspSol sequence
    tspSol.sequence.pop_back();
    tspSol.sequence.push_back(nodeInsertionInfo.insertedNode);
    tspSol.sequence.push_back(1);
}

//Removes selected node
//Assumes it's not repeated (which both 's' and 'CL' shouldn't be anyways)
void RemoveFromUnaddedNodes(vector<int>& unaddedNodes, int node)
{
    for(int i = 0; i < (int)unaddedNodes.size(); i++){
        if(unaddedNodes[i] == node){
            unaddedNodes.erase(unaddedNodes.begin() + i);
        }
    }
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

        SortAscendingByCost(insertionCost);

        double alpha = (double)rand() / RAND_MAX;
        int selected = rand() % ((int)ceil(alpha * insertionCost.size()));
        InsertIntoSolution(tspSol, insertionCost[selected]);

        RemoveFromUnaddedNodes(unaddedNodes, insertionCost[selected].insertedNode);
    }

    return tspSol;
}

/*BestImprovement funcs with differing methods to be used for RVND in LocalSearch()*/
//Tries to get the best swap possible in the solution
//As in the one that best minimizes the solution's cost
//"m" here means distMatrix
bool BestImprovementSwap(TspSolution* tspSol, double** m)
{
    //"Delta" as in the overall solution's cost
    double bestDelta = 0;
    int best_i, best_j;

    for(int i = 1; i < (int)tspSol->sequence.size()-1; i++){
        int vi = tspSol->sequence[i];
        int vi_next = tspSol->sequence[i + 1];
        int vi_prev = tspSol->sequence[i - 1];

        for(int j = i + 1; i < tspSol->sequence.size() - 1; j ++){
            int vj = tspSol->sequence[j];
            int vj_next = tspSol->sequence[i + 1];
            int vj_prev = tspSol->sequence[i - 1];

            //Calculating delta (cost) for curr possible swap
            double currDelta = -m[vi_prev][vi] - m[vi][vi_next] + m[vi_prev][vj]
                               + m[vj][vi_next] - m[vj_prev][vj] - m[vj][vj_next]
                               + m[vj_prev][vi] + m[vi][vj_next];

            //Cost comparison
            if(currDelta < bestDelta){
                bestDelta = currDelta;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap, only happens when overall solution cost lowers
    if(bestDelta < 0){
        swap(tspSol->sequence[best_i], tspSol->sequence[best_j]);
        tspSol->cost = tspSol->cost + bestDelta;
        return true;
    }

    return false;
}

//Increases the quality of current iteration's solution
//Does that by modifying the solution and evaluing the impact of each change
//Using the Random Variable Neighborhood Descent method
//Which just tests different neighborhood structures with a tad of randomness when choosing
//discarding whicever makes cost higher than currCost
void LocalSearch(TspSolution* tspSol, double** distMatrix)
{
    vector<int> NH_structures = {1, 2, 3, 4, 5};
    bool solImproved = false;

    while(!NH_structures.empty()){
        int rand_n = rand() % NH_structures.size();

        //Chose randomly
        switch(NH_structures[rand_n]){
            case 1:
                solImproved = BestImprovementSwap(tspSol, distMatrix);
                break;
            case 2:
                //solImproved = BestImprovement2Opt(tspSol, distMatrix);
                break;
            case 3:
                //solImproved = BestImprovementOrOpt(tspSol, distMatrix, 1);
                break;
            case 4:
                //solImproved = BestImprovementOrOpt(tspSol, distMatrix, 2);
                break;
            case 5:
                //solImproved = BestImprovementOrOpt(tspSol, distMatrix, 3);
                break;
        }
    }
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
            LocalSearch(&currIterSolution, data.getMatrixCost());

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
