#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

//Node matrix for burma14
double globalDistMatrix[14][14] = {
    {0, 153, 510, 706, 966, 581, 455, 70, 160, 372, 157, 567, 342, 398},
    {153, 0, 422, 664, 997, 598, 507, 197, 311, 479, 310, 581, 417, 376},
    {510, 422, 0, 289, 744, 390, 437, 491, 645, 880, 618, 374, 455, 211},
    {706, 664, 289, 0, 491, 265, 410, 664, 804, 1070, 768, 259, 499, 310},
    {966, 997, 744, 491, 0, 400, 514, 902, 990, 1261, 947, 418, 635, 636},
    {581, 598, 390, 265, 400, 0, 168, 522, 634, 910, 593, 19, 284, 239},
    {455, 507, 437, 410, 514, 168, 0, 389, 482, 757, 439, 163, 124, 232},
    {70, 197, 491, 664, 902, 522, 389, 0, 154, 406, 133, 508, 273, 355},
    {160, 311, 645, 804, 990, 634, 482, 154, 0, 276, 43, 623, 358, 498},
    {372, 479, 880, 1070, 1261, 910, 757, 406, 276, 0, 318, 898, 633, 761},
    {157, 310, 618, 768, 947, 593, 439, 133, 43, 318, 0, 582, 315, 464},
    {567, 581, 374, 259, 418, 19, 163, 508, 623, 898, 582, 0, 275, 221},
    {342, 417, 455, 499, 635, 284, 124, 273, 358, 633, 315, 275, 0, 247},
    {398, 376, 211, 310, 636, 239, 232, 355, 498, 761, 464, 221, 247, 0}
};

typedef struct {
    vector<int> sequence;
    double cost;
} TspSolution;

struct InsertionInfo {
    int insertedNode;
    int removedGraphEdge;
    double cost;
};

//Since i wanna use it as an oldschool dynamic array
double** CreateMatrix(int rows, int cols, double** matrix)
{
   //Matrix allocation
    matrix = new double*[rows];

    for(int i = 0; i < rows; ++i)
        matrix[i] = new double[cols];

    //Fillign in values
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            matrix[i][j] = globalDistMatrix[i][j];
    
    return matrix;
}

void DeleteMatrix(int rows, double** matrix)
{
    for(int i = 0; i < rows; ++i)
        delete[] matrix[i];
    delete[] matrix;
}

int BoundedRand(int min, int max)
{
    return min + rand() % (max - min + 1);
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
    int segmentMaxLength = ceil(dimension / 10.0);
    TspSolution disturbedSol;

    //Will mark the index of first and last elements of each subsequence
    //Used to make it so they don't overlap
    int subseq1Index_begin, subseq2Index_begin;
    int subseq1Index_end, subseq2Index_end;

    //Length of each subsequence and of the space among them
    int subseq2Length, inbetweenSubseqsLength; //subseq2Length uneeded

    //subseq1Index_begin = 2 + rand() % ((dimension - segmentMaxLength - 1) - 2 + 1);
    subseq1Index_begin = BoundedRand(1, dimension - segmentMaxLength - segmentMaxLength - 1);

    //subseq1Index_end = (subseq1Index_begin + 1) + rand() % ((subseq1Index_begin + segmentMaxLength - 1) - (subseq1Index_begin + 1) + 2); //Arithmetic Exception (Only happens in the 5-city (low) TSP, so it might be that)
    subseq1Index_end = BoundedRand(subseq1Index_begin + 1, subseq1Index_begin + segmentMaxLength - 1);

    //subseq2Index_begin = (subseq1Index_end + 1) + rand() % ((dimension - segmentMaxLength) - (subseq1Index_end + 1) + 2); //Arithmetic Exception may happen here as well, prob same reason
    subseq2Index_begin = BoundedRand(subseq1Index_end + 1, dimension - segmentMaxLength);
    
    //subseq2Index_end = (subseq2Index_begin + 1) + rand() % ((subseq2Index_begin + segmentMaxLength - 1) - (subseq2Index_begin + 1) + 1);
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

    return disturbedSol;
}

int main(void)
{
    double** distMatrix = nullptr;
    TspSolution tspSol, disturbedSol;

    distMatrix = CreateMatrix(14, 14, distMatrix);

    tspSol.sequence = {1, 5, 4, 8, 7, 13, 6, 12, 3, 10, 2, 11, 9, 1};
    tspSol.cost = CalculateSequenceCost(tspSol.sequence, distMatrix);

    disturbedSol = Disturbance(tspSol, distMatrix, 14);

    //Print both sequences
    cout << "Original sequence: ";
    for(int i = 0; i < (int)tspSol.sequence.size(); i++)
        cout << tspSol.sequence[i] << " ";
    
    cout << "\nCost: " << tspSol.cost << '\n';

    cout << "Disturbed sequence: ";
    for(int i = 0; i < (int)disturbedSol.sequence.size(); i++)
        cout << disturbedSol.sequence[i] << " ";

    cout << "\nCost: " << disturbedSol.cost << '\n';

    DeleteMatrix(14, distMatrix);

    return 0;
}