#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct InsertionInfo {
    int insertedNode;
    int removedGraphEdge;
    double cost;
};

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

int main(void)
{
    vector<int> unaddedNodes = {2, 3, 5, 8, 9};

    printf("Unadded Nodes (before):");
    for(int node : unaddedNodes){
        printf(" %d ", node);
    }

    RemoveFromUnaddedNodes(unaddedNodes, 9);

    printf("\n--------------------------------------\n");
    printf("Unadded Nodes (After):");

    for(int node : unaddedNodes){
        printf(" %d ", node);
    }

    return 0;
}