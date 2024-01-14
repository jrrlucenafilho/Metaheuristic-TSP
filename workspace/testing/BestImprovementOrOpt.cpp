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

int main(void)
{
    

    return 0;
}