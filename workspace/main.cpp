#include "Data.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

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

typedef struct {
    clock_t accumulatedTime = 0;
    clock_t beginTime = 0;
    clock_t endTime = 0;
} my_time_t;

//Timers, they'll hold the sum of time in clock cycles over all iterations (for avg calc)
my_time_t buildSol_time, bestImprovSwap_time, orOpt_time, orOpt2_time, orOpt3_time, two_opt_time, disturbance_time;

my_time_t* buildSol_time_ptr = &buildSol_time; 
my_time_t* bestImprovSwap_time_ptr = &bestImprovSwap_time;
my_time_t* orOpt_time_ptr = &orOpt_time;
my_time_t* orOpt2_time_ptr = &orOpt2_time;
my_time_t* orOpt3_time_ptr = &orOpt3_time;
my_time_t* two_opt_time_ptr = &two_opt_time;
my_time_t* disturbance_time_ptr = &disturbance_time;


//Prints the time for each Neighborhood Structure
void PrintNBSTimers()
{
    cout << "\tBuildSolution time: " << buildSol_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tBestImprovementSwap time: " << bestImprovSwap_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt time: " << orOpt_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt2 time: " << orOpt2_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt3 time: " << orOpt3_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\t2-Opt time: " << two_opt_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tDisturbance time: " << disturbance_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
}

//BuildSolution Utility funcs
/**
 * @brief Calcs the cost of inserting the new node into the graph. (vertex and node used interchangaebly here)
 * @param tspSol To-be-evaluated TSP Solution
 * @param unaddedVertices Vertices / Nodes that still don't belong to the solution subtour
 * @param distMatrix Matrix holding distances among graph nodes
 * @return insertionCost Cost of inserting nodes into TSP Solution graph
 **/
vector<InsertionInfo> CalcNodeInsertionCost(TspSolution& tspSol, vector<int>& unaddedVertices, double** distMatrix)
{
    vector<InsertionInfo> insertionCost = vector<InsertionInfo>((tspSol.sequence.size() - 1) * unaddedVertices.size());
    int l = 0;
    int i, j;

    //Iterating among all nodes in the solution sequence, going edge-by-edge through curr solution (pair of 2 nodes on each edge)
    for(int node1 = 0; node1 < (int)tspSol.sequence.size() - 1; node1++){
        i = tspSol.sequence[node1];
        j = tspSol.sequence[node1 + 1];

        //For each node/vertex
        for(int k : unaddedVertices){
            insertionCost[l].cost = distMatrix[i][k] + distMatrix[j][k] - distMatrix[i][j];
            insertionCost[l].insertedNode = k;
            insertionCost[l].removedGraphEdge = node1;
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
    bool nodeFound = false;

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

//Cost sorting funcs
bool CompareByCost(InsertionInfo prevInsertionInfo, InsertionInfo currInsertionInfo)
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
    tspSol.sequence.insert(tspSol.sequence.begin() + nodeInsertionInfo.removedGraphEdge + 1, nodeInsertionInfo.insertedNode);
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

double CalculateSequenceCost(vector<int>& sequence, double** m)
{
    double cost = 0;

    for(int i = 0, j = 1; i < (int)sequence.size() - 1; i++, j++){
        cost += m[sequence[i]][sequence[j]];
    }

    return cost;
}

/**
 * @brief Builds a fair solution (though still far from optimized) Using Greedy Randomized Adaptive Search (Insertion-by-cheapest)
 * @return TspSolution
 */
TspSolution BuildSolution(double** distMatrix, int dimension)
{
    buildSol_time_ptr->beginTime = std::clock();

    TspSolution tspSol;
    vector<int> addedNodes;
    tspSol.sequence = Choose3RandNodes(dimension, addedNodes);   //gets s
    vector<int> unaddedNodes = GetUnchosenNodes(dimension, addedNodes);   //gets CL

    //Calculating cost for the initial solution
    tspSol.cost = CalculateSequenceCost(tspSol.sequence, distMatrix);

    while(!unaddedNodes.empty()){
        vector<InsertionInfo> insertionCost = CalcNodeInsertionCost(tspSol, unaddedNodes, distMatrix);

        SortAscendingByCost(insertionCost);

        double alpha = ((double)rand() + 1) / RAND_MAX; //rand() + 1 so that alpha can't be zero, which throws an Arithmetic exception on next line
        int selected = rand() % ((int)ceil(alpha * insertionCost.size()));
        InsertIntoSolution(tspSol, insertionCost[selected]);

        RemoveFromUnaddedNodes(unaddedNodes, insertionCost[selected].insertedNode);
    }

    tspSol.cost = CalculateSequenceCost(tspSol.sequence, distMatrix);

    buildSol_time_ptr->endTime = std::clock();
    buildSol_time_ptr->accumulatedTime += buildSol_time_ptr->endTime - buildSol_time_ptr->beginTime;

    return tspSol;
}

/*BestImprovement funcs with differing methods to be used for RVND in LocalSearch()*/
//Tries to get the best swap possible in the solution
//As in the one that best minimizes the solution's cost
//"m" here means distMatrix
bool BestImprovementSwap(TspSolution& tspSol, double** m, int dimension)
{
    bestImprovSwap_time_ptr->beginTime = std::clock();

    //"Delta" as in the overall solution's cost
    double bestDelta = 0;
    double currDelta;
    double toBeRemovedSection;
    int best_i = 0, best_j = 0;
    int graphSize = dimension + 1;

    for(int i = 1; i < graphSize - 1; i++){
        for(int j = i + 1; j < graphSize - 1; j++){
            toBeRemovedSection = -(m[tspSol.sequence[i]][tspSol.sequence[i + 1]] + m[tspSol.sequence[i]][tspSol.sequence[i - 1]] 
                                 + m[tspSol.sequence[j]][tspSol.sequence[j + 1]] + m[tspSol.sequence[j]][tspSol.sequence[j - 1]]);
        
            //In case of neighboring nodes
            if(j == i + 1){
                currDelta = - m[tspSol.sequence[i - 1]][tspSol.sequence[i]] - m[tspSol.sequence[j]][tspSol.sequence[j + 1]] 
                            + m[tspSol.sequence[i - 1]][tspSol.sequence[j]] + m[tspSol.sequence[i]][tspSol.sequence[j + 1]];
            }else{
                currDelta = toBeRemovedSection + m[tspSol.sequence[i]][tspSol.sequence[j + 1]] + m[tspSol.sequence[i]][tspSol.sequence[j - 1]] 
                            + m[tspSol.sequence[j]][tspSol.sequence[i + 1]] + m[tspSol.sequence[j]][tspSol.sequence[i - 1]];
            }

            if(currDelta < bestDelta){
                bestDelta = currDelta;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap, only happens when overall solution cost lowers
    if(bestDelta < 0){
        swap(tspSol.sequence[best_i], tspSol.sequence[best_j]);
        tspSol.cost = tspSol.cost + bestDelta;
        return true;

        bestImprovSwap_time_ptr->endTime = std::clock();
        bestImprovSwap_time_ptr->accumulatedTime += bestImprovSwap_time_ptr->endTime - bestImprovSwap_time_ptr->beginTime;
    }

    bestImprovSwap_time_ptr->endTime = std::clock();
    bestImprovSwap_time_ptr->accumulatedTime += bestImprovSwap_time_ptr->endTime - bestImprovSwap_time_ptr->beginTime;

    return false;
}

bool BestImprovement2Opt(TspSolution& tspSol, double** m, int dimension)
{
    two_opt_time_ptr->beginTime = std::clock();

    double bestDelta = 0;
    double initialDelta, currDelta;
    int best_i = 0, best_j = 0;
    int graphSize = dimension + 1;

    for(int i = 1; i < graphSize - 1; i++){
        initialDelta = -m[tspSol.sequence[i - 1]][tspSol.sequence[i]];

        for(int j = i + 1; j < graphSize - 1; j++){
            currDelta = initialDelta -m[tspSol.sequence[j]][tspSol.sequence[j + 1]] 
                                     + m[tspSol.sequence[i - 1]][tspSol.sequence[j]] 
                                     + m[tspSol.sequence[i]][tspSol.sequence[j + 1]];

            if(currDelta < bestDelta){
                bestDelta = currDelta;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap
    if(bestDelta < 0){
        reverse(tspSol.sequence.begin() + best_i, tspSol.sequence.begin() + best_j + 1);
        tspSol.cost = tspSol.cost + bestDelta;

        two_opt_time_ptr->endTime = std::clock();
        two_opt_time_ptr->accumulatedTime += two_opt_time_ptr->endTime - two_opt_time_ptr->beginTime;

        return true;
    }

    two_opt_time_ptr->endTime = std::clock();
    two_opt_time_ptr->accumulatedTime += two_opt_time_ptr->endTime - two_opt_time_ptr->beginTime;

    return false;
}

bool BestImprovementOrOpt(TspSolution& tspSol, double** m, int dimension, int movedBlockSize)
{
    //Just to account for var creation time, prob quite smallish
    my_time_t var_creation_time;
    var_creation_time.beginTime = std::clock();

    double bestDelta = 0;
    double initialDelta, currDelta = 0;
    int best_i = 0, best_j = 0;
    int graphSize = dimension + 1;

    var_creation_time.endTime = std::clock();
    var_creation_time.accumulatedTime += var_creation_time.endTime - var_creation_time.beginTime;

    //Reinsertion Case
    if(movedBlockSize == 1){
        orOpt_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 1; i++){
            initialDelta = -m[tspSol.sequence[i - 1]][tspSol.sequence[i]] 
                           - m[tspSol.sequence[i]][tspSol.sequence[i + 1]] 
                           + m[tspSol.sequence[i - 1]][tspSol.sequence[i + 1]];
        
            for(int j = 1; j < graphSize - 1; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol.sequence[j]][tspSol.sequence[j + 1]]
                                                 + m[tspSol.sequence[i]][tspSol.sequence[j]] 
                                                 + m[tspSol.sequence[i]][tspSol.sequence[j + 1]];
                    }else{
                        currDelta = initialDelta -m[tspSol.sequence[j]][tspSol.sequence[j - 1]] 
                                                 + m[tspSol.sequence[i]][tspSol.sequence[j]] 
                                                 + m[tspSol.sequence[j - 1]][tspSol.sequence[i]];
                    }

                    if(currDelta < bestDelta){
                        bestDelta = currDelta;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        orOpt_time_ptr->endTime = std::clock();
        orOpt_time_ptr->accumulatedTime += orOpt_time_ptr->endTime - orOpt_time_ptr->beginTime;
    }

    //OrOpt-2 case
    if(movedBlockSize == 2){
        orOpt2_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 2; i++){
            initialDelta = -m[tspSol.sequence[i - 1]][tspSol.sequence[i]]
                           - m[tspSol.sequence[i + 1]][tspSol.sequence[i + 2]] 
                           + m[tspSol.sequence[i - 1]][tspSol.sequence[i + 2]];

            for(int j = 1; j < graphSize - 3; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol.sequence[j + 1]][tspSol.sequence[j + 2]] 
                                                 + m[tspSol.sequence[i + 1]][tspSol.sequence[j + 2]] 
                                                 + m[tspSol.sequence[i]][tspSol.sequence[j + 1]];
                    }else{
                        currDelta = initialDelta -m[tspSol.sequence[j - 1]][tspSol.sequence[j]] 
                                                 + m[tspSol.sequence[j - 1]][tspSol.sequence[i]] 
                                                 + m[tspSol.sequence[i + 1]][tspSol.sequence[j]];
                    }

                    if(currDelta < bestDelta){
                        bestDelta = currDelta;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        orOpt2_time_ptr->endTime = std::clock();
        orOpt2_time_ptr->accumulatedTime += orOpt2_time_ptr->endTime - orOpt2_time_ptr->beginTime;
    }

    //OrOpt-3 case
    if(movedBlockSize == 3){
        orOpt3_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 3; i++){
            initialDelta = -m[tspSol.sequence[i - 1]][tspSol.sequence[i]]
                           - m[tspSol.sequence[i + 2]][tspSol.sequence[i + 3]] 
                           + m[tspSol.sequence[i - 1]][tspSol.sequence[i + 3]];
    
            for(int j = 1; j < graphSize - 4; j++){
                if(i != j){
                    if(i < j){
                        currDelta = initialDelta -m[tspSol.sequence[j + 2]][tspSol.sequence[j + 3]] 
                                                 + m[tspSol.sequence[i]][tspSol.sequence[j + 2]]
                                                 + m[tspSol.sequence[i + 2]][tspSol.sequence[j + 3]];
                    }else{
                        currDelta = initialDelta -m[tspSol.sequence[j - 1]][tspSol.sequence[j]] 
                                                 + m[tspSol.sequence[j - 1]][tspSol.sequence[i]] 
                                                 + m[tspSol.sequence[i + 2]][tspSol.sequence[j]];
                    }

                    if(currDelta < bestDelta){
                        bestDelta = currDelta;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        orOpt3_time_ptr->endTime = std::clock();
        orOpt3_time_ptr->accumulatedTime += orOpt3_time_ptr->endTime - orOpt3_time_ptr->beginTime;
    }

    if(bestDelta < 0){
        //Reinsertion Case
        if(movedBlockSize == 1){
            orOpt_time_ptr->beginTime = std::clock();

            int reinsertedNode = tspSol.sequence[best_i];
            tspSol.sequence.erase(tspSol.sequence.begin() + best_i);
            tspSol.sequence.insert(tspSol.sequence.begin() + best_j, reinsertedNode);

            tspSol.cost += bestDelta;

            orOpt_time_ptr->endTime = std::clock();
            orOpt_time_ptr->accumulatedTime += orOpt_time_ptr->endTime - orOpt_time_ptr->beginTime;
            orOpt_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;
        }

        //OrOpt-2 case
        if(movedBlockSize == 2){
            orOpt2_time_ptr->beginTime = std::clock();

            vector<int> reinsertSequence(tspSol.sequence.begin() + best_i, tspSol.sequence.begin() + best_i + 2);
            tspSol.sequence.erase(tspSol.sequence.begin() + best_i, tspSol.sequence.begin() + best_i + 2);
            tspSol.sequence.insert(tspSol.sequence.begin() + best_j, reinsertSequence.begin(), reinsertSequence.end());

            tspSol.cost += bestDelta;

            orOpt2_time_ptr->endTime = std::clock();
            orOpt2_time_ptr->accumulatedTime += orOpt2_time_ptr->endTime - orOpt2_time_ptr->beginTime;
            orOpt2_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;
        }

        //OrOpt-3 case
        if(movedBlockSize == 3){
            orOpt3_time_ptr->beginTime = std::clock();

            vector<int> reinsertSequence(tspSol.sequence.begin() + best_i, tspSol.sequence.begin() + best_i + 3);
            tspSol.sequence.erase(tspSol.sequence.begin() + best_i, tspSol.sequence.begin() + best_i + 3);
            tspSol.sequence.insert(tspSol.sequence.begin() + best_j, reinsertSequence.begin(), reinsertSequence.end());

            tspSol.cost += bestDelta;

            orOpt3_time_ptr->endTime = std::clock();
            orOpt3_time_ptr->accumulatedTime += orOpt3_time_ptr->endTime - orOpt3_time_ptr->beginTime;
            orOpt3_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;
        }
        
        return true;
    }

    return false;
}

//Increases the quality of current iteration's solution
//Does that by modifying the solution and evaluing the impact of each change
//Using the Random Variable Neighborhood Descent method
//Which just tests different neighborhood structures with a tad of randomness when choosing
//discarding whichever makes cost higher than currCost
void LocalSearch(TspSolution& tspSol, double** distMatrix, int dimension)
{
    vector<int> NH_structures = {1, 2, 3, 4, 5};
    bool solutionImproved = false;

    while(!NH_structures.empty()){
        int rand_n = rand() % NH_structures.size();

        //Chooses randomly
        switch(NH_structures[rand_n]){
            case 1:
                solutionImproved = BestImprovementSwap(tspSol, distMatrix, dimension);
                break;
            case 2:
                solutionImproved = BestImprovement2Opt(tspSol, distMatrix, dimension);
                break;
            case 3:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, dimension, 1);
                break;
            case 4:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, dimension, 2);
                break;
            case 5:
                solutionImproved = BestImprovementOrOpt(tspSol, distMatrix, dimension, 3);
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

int BoundedRand(int min, int max)
{
    return min + rand() % (max - min + 1);
}

//Return TspSolution, may have to call a cost-calc function inside here
TspSolution Disturbance(TspSolution& tspSol, double** m, int dimension)
{
    disturbance_time_ptr->beginTime = std::clock();

    vector<int> copiedSequence = tspSol.sequence;
    int segmentMaxLength = ceil(dimension / 10.0);
    TspSolution disturbedSol;

    //Will mark the index of first and last elements of each subsequence
    //Used to make it so they don't overlap
    int subseq1Index_begin, subseq2Index_begin;
    int subseq1Index_end, subseq2Index_end;

    //Length of each subsequence and of the space among them
    int subseq2Length, inbetweenSubseqsLength; //subseq2Length uneeded

    subseq1Index_begin = BoundedRand(1, dimension - segmentMaxLength - segmentMaxLength - 1);
    subseq1Index_end = BoundedRand(subseq1Index_begin + 1, subseq1Index_begin + segmentMaxLength - 1);
    subseq2Index_begin = BoundedRand(subseq1Index_end + 1, dimension - segmentMaxLength);
    subseq2Index_end = BoundedRand(subseq2Index_begin + 1, subseq2Index_begin + segmentMaxLength - 1);

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

    //Calculating cost and attributing to disturbed solution
    disturbedSol.cost = CalculateSequenceCost(copiedSequence, m);
    disturbedSol.sequence = copiedSequence;

    disturbance_time_ptr->endTime = std::clock();
    disturbance_time_ptr->accumulatedTime += disturbance_time_ptr->endTime - disturbance_time_ptr->beginTime;

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
            //Tries to enhance the fairly-guessed solution
            //By doing small modifications to it
            LocalSearch(currIterSolution, data.getMatrixCost(), data.getDimension());

            if(currIterSolution.cost < currBestSolution.cost){
                currBestSolution = currIterSolution;
                iterILS = 0;
            }

            //If not possible to make it better, shake the current best solution up a lil'
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
    double costsSum = 0;
    auto data = Data(argc, argv[1]);
    TspSolution tspSol;
    my_time_t ILS_iter_time;

    data.read();
    data.reformatMatrix();
    size_t n = data.getDimension();

    cout << "Dimension: " << n << '\n';
    //cout << "DistanceMatrix: " << '\n';
    //data.printMatrixDist();
    cout << "Wait for it...\n";

    //Defining Iters
    if(data.getDimension() >= 150){
        maxIterILS = data.getDimension() / 2.0;
    }else{
        maxIterILS = data.getDimension();
    }

    cout << "-------------------------------\n";
    cout << "10 Iteration costs:\n";

    //10 execs for the sum avg calc
    for(int i = 0; i < 10; i++){
        srand(static_cast<unsigned int>(time(0)));

        ILS_iter_time.beginTime = std::clock();

        tspSol = IteratedLocalSearch(maxIter, maxIterILS, data);

        tspSol.cost = CalculateSequenceCost(tspSol.sequence, data.getMatrixCost());
        costsSum += tspSol.cost;

        ILS_iter_time.endTime = std::clock();
        ILS_iter_time.accumulatedTime += ILS_iter_time.endTime - ILS_iter_time.beginTime;

        cout << "Cost of s (iter " << i + 1 << "): " << tspSol.cost << '\n';
    }

    cout << "-------------------------------\n";
    cout << "Instance Name: " << data.getInstanceName() << '\n';
    cout << "Solution s = ";

    for(size_t i = 0; i < tspSol.sequence.size() - 1; i++){
        cout << tspSol.sequence[i] << " -> ";
    }
    cout << "1\n";
    cout << "Average s cost: " << costsSum / 10 << '\n';
    cout << "Average CPU execution time: " << (ILS_iter_time.accumulatedTime / static_cast<double>(CLOCKS_PER_SEC)) / 10 << " s \n";
    cout << "-------------------------------\n";
    cout << "Average Times on each Neighborhood Structure:\n";
    PrintNBSTimers();

    return 0;
}
