#include "Data.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[])
{
    Data data = Data(argc, argv[1]);
    data.read();

    data.printMatrixDist();

    data.reformatMatrix();

    cout << "\n---------Reformatted matrix---------\n" << '\n';
    cout << "Dimension: " << data.getDimension() << '\n';

    data.printMatrixDist();
    return 0;
}