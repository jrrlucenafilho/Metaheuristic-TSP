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

//Chooses 3 random nodes in the graph from 1 up to dimension (numOfNodes)
//Makes a sequence from them following {1, x, y, z, 1} and puts the 3 added nodes into solAddedNodes
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

//Inserts a node into the tsp solution
void InsertIntoSolution(TspSolution& tspSol, InsertionInfo& nodeInsertionInfo)
{
    //First adds up the cost of the to-be-added node
    tspSol.cost += nodeInsertionInfo.cost;

    //Added the node into the tspSol sequence
    tspSol.sequence.pop_back();
    tspSol.sequence.push_back(nodeInsertionInfo.insertedNode);
    tspSol.sequence.push_back(1);
}

int main(void)
{
    vector<int> starter_seq;
    vector<int> solAddedNodes;

    starter_seq = Choose3RandNodes(10, solAddedNodes);

    printf("Solution Nodes:");
    for(int i = 0; i < (int)starter_seq.size(); i++){
        printf(" %d ", starter_seq[i]);
    }

    //Defining a tspSol
    TspSolution tspSol;
    tspSol.sequence = starter_seq;
    tspSol.cost = 1;

    //Defining an insertionInfo
    InsertionInfo insInfo;
    insInfo.cost = 20;
    insInfo.insertedNode = 3;
    insInfo.removedGraphEdge = 3;

    InsertIntoSolution(tspSol, insInfo);

    printf("\n-------------------------------------------\n");
    printf("Pós-Inserção: ");

    for(int i = 0; i < (int)tspSol.sequence.size(); i++){
        printf(" %d ", tspSol.sequence[i]);
    }

    printf("\nCost: %lf\n", tspSol.cost);

    return 0;
}