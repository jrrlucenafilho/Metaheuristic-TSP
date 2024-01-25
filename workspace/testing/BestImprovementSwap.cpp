#include <iostream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef struct {
    vector<int> sequence;
    double cost;
} TspSolution;

bool BestImprovementSwap(TspSolution* tspSol, double** m, int dimension)
{
    //"Delta" as in the overall solution's cost
    double bestDelta = 0;
    int best_i = 0, best_j = 0;

    //Iterates over all pairs of vertices in the tour. For each pair (i, j), it calculates the change in cost (currDelta) that would result from swapping these two vertices.
    for(int i = 1; i < dimension - 1; i++){
        int vi = tspSol->sequence[i];
        int vi_next = tspSol->sequence[i + 1];
        int vi_prev = tspSol->sequence[i - 1];

        for(int j = i + 1; j < dimension - 1; j++){
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

int main(void)
{
    
}