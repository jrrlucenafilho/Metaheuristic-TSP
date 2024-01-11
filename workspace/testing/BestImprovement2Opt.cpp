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
    double cost;
};

//2-Opt (TOTEST)
bool BestImprovement2Opt(TspSolution* tspSol, double** m)
{
    //"Delta" as in the overall solution's cost
    double bestDelta = 0;
    double initialDelta, currDelta;
    int best_i = 0, best_j = 0;

    for(int i = 1; i < (int)tspSol->sequence.size() - 1; i++){
        initialDelta = -(m[tspSol->sequence[i-1]][tspSol->sequence[i]]);

        for(int j = i + 1; j < (int)tspSol->sequence.size() - 1; j++){
            currDelta = initialDelta - m[tspSol->sequence[j]][tspSol->sequence[j + 1]] 
            + m[tspSol->sequence[i - 1]][tspSol->sequence[j]] 
            + m[tspSol->sequence[i]][tspSol->sequence[j + 1]];

            //Cost comparison
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

int main(void)
{
    

    return 0;
}