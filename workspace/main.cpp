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

vector<int> Choose3RandNodes(double** distMatrix)
{
    //TODO: What the func says duh. Test it too
    //Gotta make the next func know which ones haven't been chosen
    //Gotta have size of graph

    vector<int> starterSequence;
    int randChosenEdge;     //Randomly chooses an edge in 2d array
    bool randChosenNode;    //true -> chosses up node, false does opposite

    starterSequence.push_back(1);

    //Algo for now is: Choose randomly from (0 to numOfEdges-1)
    //And then randomly assign true or false for randChosenNode, to +1 or -1 on chosen edge
}

/**
 * @brief Builds a fair solution (though still far from optimized) Using Greedy Randomized Adaptive Search (Insertion-by-cheapest)
 * @return TspSolution
 */
TspSolution BuildSolution(double** distMatrix)
{
    TspSolution tspSol;
    tspSol.sequence = Choose3RandNodes(distMatrix);   //gets s
    //vector<int> unaddedNodes = getUnchosenNodes();   //gets CL

    while(!unaddedNodes.empty()){
        vector<InsertionInfo> insertionCost = CalcNodeInsertionCost(tspSol, unaddedNodes, distMatrix);

        insertionCost;

        double alpha = (double)rand() / RAND_MAX;

        int selected = rand() % ((int)ceil(alpha * insertionCost.size()));

        //insertIntoSolution(tspSol, insertionCost[selected].insertedNode); //k
    }

    return tspSol;
}

TspSolution IteratedLocalSearch(int maxIters, int maxIterILS, double** distMatrix)
{
    TspSolution bestOfAllSolution;
    bestOfAllSolution.cost = INFINITY;

    for(int i = 0; i < maxIters; i++){
        //Builds a beginning solution based on fair guesses
        TspSolution currIterSolution = BuildSolution(distMatrix);
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
    IteratedLocalSearch(5, 5, data.getMatrixCost());

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
