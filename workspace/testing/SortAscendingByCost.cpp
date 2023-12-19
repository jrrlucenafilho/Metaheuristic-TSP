#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct InsertionInfo {
    int insertedNode;
    int removedGraphEdge;
    double cost;
};

//Sorting the costs
bool CompareByCost(InsertionInfo& prevInsertionInfo, InsertionInfo& currInsertionInfo)
{
    return prevInsertionInfo.cost < currInsertionInfo.cost;
}

void SortAscendingByCost(vector<InsertionInfo>& insertionInfo)
{
    sort(insertionInfo.begin(), insertionInfo.end(), CompareByCost);
}

int main(void)
{
    vector<InsertionInfo> insertionInfo;
    
    //Creating InsertionInfos
    InsertionInfo if1;
    InsertionInfo if2;
    InsertionInfo if3;
    InsertionInfo if4;
    InsertionInfo if5;
    InsertionInfo if6;
    InsertionInfo if7;
    InsertionInfo if8;
    InsertionInfo if9;
    InsertionInfo if10;

    if1.cost = 23;
    if2.cost = 30;
    if3.cost = 5;
    if4.cost = 99;
    if5.cost = 2300;
    if6.cost = 49;
    if7.cost = 300;
    if8.cost = 230;
    if9.cost = 100;
    if10.cost = 10000;

    insertionInfo.push_back(if1);
    insertionInfo.push_back(if2);
    insertionInfo.push_back(if3);
    insertionInfo.push_back(if4);
    insertionInfo.push_back(if5);
    insertionInfo.push_back(if6);
    insertionInfo.push_back(if7);
    insertionInfo.push_back(if8);
    insertionInfo.push_back(if9);
    insertionInfo.push_back(if10);

    SortAscendingByCost(insertionInfo);

    printf("Ordered costs: \n");
    for(InsertionInfo k : insertionInfo){
        printf(" %lf ", k.cost);
    }

    return 0;
}