//NOTE: may have to change all tspSol->sequence.size() to "dimension - 1" on the neighborhood structure functions

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

        //Handles if rand chooses 0 or 1, since distances 1-indexed here and 1 is the first and last city
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
            break;
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
    int best_i = 0, best_j = 0;

    //Iterates over all pairs of vertices in the tour. For each pair (i, j), it calculates the change in cost (currDelta) that would result from swapping these two vertices.
    for(int i = 1; i < (int)tspSol->sequence.size() - 1; i++){
        int vi = tspSol->sequence[i];
        int vi_next = tspSol->sequence[i + 1];
        int vi_prev = tspSol->sequence[i - 1];

        for(int j = i + 1; j < (int)tspSol->sequence.size() - 1; j++){
            int vj = tspSol->sequence[j];
            int vj_next = tspSol->sequence[j + 1];
            int vj_prev = tspSol->sequence[j - 1];

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

//2-Opt (TOTEST)
bool BestImprovement2Opt(TspSolution* tspSol, double** m)
{
    double bestDelta = 0;
    double initialDelta, currDelta;
    int best_i = 0, best_j = 0;

    for(int i = 1; i < (int)tspSol->sequence.size() - 1; i++){
        initialDelta = -m[tspSol->sequence[i - 1]][tspSol->sequence[i]];

        for(int j = i + 1; j < (int)tspSol->sequence.size() - 1; j++){
            currDelta = initialDelta -m[tspSol->sequence[j]][tspSol->sequence[j + 1]] 
                                     + m[tspSol->sequence[i - 1]][tspSol->sequence[j]] 
                                     + m[tspSol->sequence[i]][tspSol->sequence[j + 1]];

            if(currDelta < bestDelta){
                bestDelta = currDelta;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap
    if(bestDelta < 0){
        reverse(tspSol->sequence.begin() + best_i, tspSol->sequence.begin() + best_j + 1);
        tspSol->cost = tspSol->cost + bestDelta;
        return true;
    }

    return false;
}

//OrOpt (TOTEST)
bool BestImprovementOrOpt(TspSolution* tspSol, double** m, int movedBlockSize)
{
    double bestDelta = 0;
    double initialDelta, currDelta;
    int best_i = 0, best_j = 0;

    //Reinsertion Case
    if(movedBlockSize == 1){
        for(int i = 1; i < (int)tspSol->sequence.size(); i++){
            initialDelta = -m[tspSol->sequence[i - 1]][tspSol->sequence[i]] 
                           - m[tspSol->sequence[i]][tspSol->sequence[i + 1]] 
                           + m[tspSol->sequence[i - 1]][tspSol->sequence[i + 1]];
        
            for(int j = 1; j < (int)tspSol->sequence.size() - 1; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol->sequence[j]][tspSol->sequence[j + 1]] 
                                                 + m[tspSol->sequence[i]][tspSol->sequence[j]] 
                                                 + m[tspSol->sequence[i]][tspSol->sequence[j + 1]];
                    }else{
                        currDelta = initialDelta -m[tspSol->sequence[j]][tspSol->sequence[j - 1]] 
                                                 + m[tspSol->sequence[i]][tspSol->sequence[j]] 
                                                 + m[tspSol->sequence[j - 1]][tspSol->sequence[i]];
                    }
                }

                if(currDelta < bestDelta){
                    bestDelta = currDelta;
                    best_i = i;
                    best_j = j;
                }
            }
        }
    }

    //OrOpt-2 case
    if(movedBlockSize == 2){
        for(int i = 1; i < (int)tspSol->sequence.size() - 2; i++){
            initialDelta = -m[tspSol->sequence[i - 1]][tspSol->sequence[i]]
                           - m[tspSol->sequence[i + 1]][tspSol->sequence[i + 2]] 
                           + m[tspSol->sequence[i - 1]][tspSol->sequence[i + 2]];

            for(int j = 1; j < (int)tspSol->sequence.size() - 3; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol->sequence[j + 1]][tspSol->sequence[j + 2]] 
                                                 + m[tspSol->sequence[i + 1]][tspSol->sequence[j + 2]] 
                                                 + m[tspSol->sequence[i]][tspSol->sequence[j + 1]];
                    }else{
                        currDelta = initialDelta -m[tspSol->sequence[j - 1]][tspSol->sequence[j]] 
                                                 + m[tspSol->sequence[j - 1]][tspSol->sequence[i]] 
                                                 + m[tspSol->sequence[i + 1]][tspSol->sequence[j]];
                    }

                    if(currDelta < bestDelta){
                        bestDelta = currDelta;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
    }

    //OrOpt-3 case
    if(movedBlockSize == 3){
        for(int i = 1; i < (int)tspSol->sequence.size() - 3; i++){
            initialDelta = -m[tspSol->sequence[i - 1]][tspSol->sequence[i]]
                           - m[tspSol->sequence[i + 2]][tspSol->sequence[i + 3]] 
                           + m[tspSol->sequence[i - 1]][tspSol->sequence[i + 3]];
    
            for(int j = 1; j < (int)tspSol->sequence.size() - 4; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol->sequence[j + 2]][tspSol->sequence[j + 3]] 
                                                 + m[tspSol->sequence[i]][tspSol->sequence[j + 2]]
                                                 + m[tspSol->sequence[i + 2]][tspSol->sequence[j + 3]];
                    }else{
                        currDelta = initialDelta -m[tspSol->sequence[j - 1]][tspSol->sequence[j]] 
                                                 + m[tspSol->sequence[j - 1]][tspSol->sequence[i]] 
                                                 + m[tspSol->sequence[i + 2]][tspSol->sequence[j]];
                    }
                }

                if(currDelta < bestDelta){
                    bestDelta = currDelta;
                    best_i = i;
                    best_j = j;
                }
            }
        }
    }

    if(bestDelta < 0){
        //Reinsertion Case
        if(movedBlockSize == 1){
            int reinsertedNode = tspSol->sequence[best_i];
            tspSol->sequence.erase(tspSol->sequence.begin() + best_i);
            tspSol->sequence.insert(tspSol->sequence.begin() + best_j, reinsertedNode);
        }

        //OrOpt-2 case
        if(movedBlockSize == 2){
            vector<int> reinsertSequence(tspSol->sequence.begin() + best_i, tspSol->sequence.begin() + best_i + 2);
            tspSol->sequence.erase(tspSol->sequence.begin() + best_i, tspSol->sequence.begin() + best_i + 2);
            tspSol->sequence.insert(tspSol->sequence.begin() + best_j, reinsertSequence.begin(), reinsertSequence.end());
        }

        //OrOpt-3 case
        if(movedBlockSize == 3){
            vector<int> reinsertSequence(tspSol->sequence.begin() + best_i, tspSol->sequence.begin() + best_i + 3);
            tspSol->sequence.erase(tspSol->sequence.begin() + best_i, tspSol->sequence.begin() + best_i + 3);
            tspSol->sequence.insert(tspSol->sequence.begin() + best_j, reinsertSequence.begin(), reinsertSequence.end());
        }

        tspSol->cost += bestDelta;
        
        return true;
    }

    return false;
}

//Increases the quality of current iteration's solution
//Does that by modifying the solution and evaluing the impact of each change
//Using the Random Variable Neighborhood Descent method
//Which just tests different neighborhood structures with a tad of randomness when choosing
//discarding whichever makes cost higher than currCost
void LocalSearch(TspSolution* tspSol, double** distMatrix)
{
    vector<int> NH_structures = {1, 2, 3, 4, 5};
    bool solutionImproved = false;

    while(!NH_structures.empty()){
        int rand_n = rand() % NH_structures.size();

        //Chooses randomly
        switch(NH_structures[rand_n]){
            case 1:
                solutionImproved = BestImprovementSwap(tspSol, distMatrix);
                break;
            case 2:
                solutionImproved = BestImprovement2Opt(tspSol, distMatrix);
                break;
            case 3:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, 1);
                break;
            case 4:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, 2);
                break;
            case 5:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, 3);
                break;
        }

        //Checks if solution's improved, erasing neighborhood structure if not
        if(solutionImproved){
            NH_structures = {1, 2, 3, 4, 5};
        }else{
            NH_structures.erase(NH_structures.begin() + rand_n);
        }
    }
}

double CalculateSequenceCost(vector<int>& sequence, double** m)
{
    double cost = 0;

    for(int i = 0, j = 1; i < (int)sequence.size() - 1; i++, j++){
        cost += m[sequence[i]][sequence[j]];
    }

    return cost;
}

//Return TspSolution, may have to call a cost-calc function inside here
TspSolution Disturbance(TspSolution& tspSol, double** m, int dimension)
{
    vector<int> copiedSequence = tspSol.sequence;
    int segmentMaxLength = dimension / 10;
    TspSolution disturbedSol;

    //Will mark the index of first and last elements of each subsequence
    //Used to make it so they don't overlap
    int subseq1Index_begin, subseq2Index_begin;
    int subseq1Index_end, subseq2Index_end;

    //Length of each subsequence and of the space among them
    int subseq2Length, inbetweenSubseqsLength; //subseq2Length uneeded

    subseq1Index_begin = 2 + rand() % ((dimension - segmentMaxLength - 1) - 2 + 1);
    subseq1Index_end = (subseq1Index_begin + 1) + rand() % ((subseq1Index_begin + segmentMaxLength - 1) - (subseq1Index_begin + 1) + 1);

    subseq2Index_begin = (subseq1Index_end + 1) + rand() % ((dimension - segmentMaxLength) - (subseq1Index_end + 1) + 1);
    subseq2Index_end = (subseq2Index_begin + 1) + rand() % ((subseq2Index_begin + segmentMaxLength - 1) - (subseq2Index_begin + 1) + 1);
    
    //Actually making the subsequences from rand calc'd indexes
    vector<int> subseq1(copiedSequence.begin() + subseq1Index_begin, copiedSequence.begin() + subseq1Index_end);
    vector<int> subseq2(copiedSequence.begin() + subseq2Index_begin, copiedSequence.begin() + subseq2Index_end);

    //Lengths and space calc
    //subseq1Length = subseq1Index_end - subseq1Index_begin;
    subseq2Length = subseq2Index_end - subseq2Index_begin;
    inbetweenSubseqsLength = subseq2Index_begin - subseq1Index_end;

    //Disturbing copied sequence
    //Swapping these not-necessarily-equally-sized subseqs in the main sequence
    copiedSequence.erase(copiedSequence.begin() + subseq1Index_begin, copiedSequence.begin() + subseq1Index_end);
    copiedSequence.insert(copiedSequence.begin() + subseq1Index_begin, subseq2.begin(), subseq2.end());

    copiedSequence.erase(copiedSequence.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength, copiedSequence.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength + subseq2Length);
    copiedSequence.insert(copiedSequence.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength, subseq1.begin(), subseq1.end());

    //Calculating and attributing to disturbed solution
    disturbedSol.cost = CalculateSequenceCost(copiedSequence, m);
    disturbedSol.sequence = copiedSequence;

    return disturbedSol;
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
            currIterSolution = Disturbance(currBestSolution, data.getMatrixCost(), data.getDimension());
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
    int maxIter = 50;
    int maxIterILS;
    auto data = Data(argc, argv[1]);

    data.read();

    size_t n = data.getDimension();

    cout << "Dimension: " << n << endl;
    cout << "DistanceMatrix: " << endl;
    data.printMatrixDist();

    //Defining Iters
    if(data.getDimension() >= 150){
        maxIterILS = data.getDimension() / 2;
    }else{
        maxIterILS = data.getDimension();
    }

    IteratedLocalSearch(maxIter, maxIterILS, data);

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
