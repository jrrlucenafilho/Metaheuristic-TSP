#include <iostream>
#include <vector>

using namespace std;

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

int main(void)
{
    vector<int> starter_seq;
    vector<int> solAddedNodes;
    vector<int> unaddedNodes;

    starter_seq = Choose3RandNodes(10, solAddedNodes);

    printf("Added Nodes:\n");
    for(int node : solAddedNodes){
        printf("Node added: %d\n", node);
    }

    printf("\ns:\n");
    for(int i = 0; i < (int)starter_seq.size(); i++){
        printf("Solution node %d: %d\n", i, starter_seq[i]);
    }

    printf("\n--------------------------------------------------\n");

    unaddedNodes = GetUnchosenNodes(10, solAddedNodes);

    printf("CL:\n");
    for(int node : unaddedNodes){
        printf("Unadded Node: %d\n", node);
    }

    return 0;
}