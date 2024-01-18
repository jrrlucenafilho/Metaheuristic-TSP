#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef struct {
    vector<int> sequence;
    double cost;
} TspSolution;

struct InsertionInfo {
    int insertedNode;
    int removedGraphEdge;
    double cost;
};

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
    int subseq1Length, subseq2Length, inbetweenSubseqsLength;

    subseq1Index_begin = 2 + rand() % ((dimension - segmentMaxLength - 1) - 2 + 1);
    subseq1Index_end = (subseq1Index_begin + 1) + rand() % ((subseq1Index_begin + segmentMaxLength - 1) - (subseq1Index_begin + 1) + 1);

    subseq2Index_begin = (subseq1Index_end + 1) + rand() % ((dimension - segmentMaxLength) - (subseq1Index_end + 1) + 1);
    subseq2Index_end = (subseq2Index_begin + 1) + rand() % ((subseq2Index_begin + segmentMaxLength - 1) - (subseq2Index_begin + 1) + 1);
    
    //Actually making the subsequences from rand calc'd indexes
    vector<int> subseq1(copiedSequence.begin() + subseq1Index_begin, copiedSequence.begin() + subseq1Index_end);
    vector<int> subseq2(copiedSequence.begin() + subseq2Index_begin, copiedSequence.begin() + subseq2Index_end);

    //Lengths and space calc
    subseq1Length = subseq1Index_end - subseq1Index_begin;
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

int main(void)
{
    

    return 0;
}