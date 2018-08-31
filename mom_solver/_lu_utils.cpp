#include<iostream>
#include"./_lu_utils.h"

using namespace boost::numeric::ublas;
using namespace std;
int main()
{
    matrix_defined b;
    std::cout << "matrix0 = " << std::endl;
    std::cout << *(b.print()) << std::endl;
    //RankRevealingMethod_tradition test;
    Result_RankRevealingMethod_tradition result;
    result = RankRevealingMethod_tradition().LUDec(b, b.size1(),b.size2(),0.01f,-1);
    std::cout << *(result.L) << std::endl;
    std::cout << *(result.U) << std::endl;
    std::cout << "matrix1 = \n";
    std::cout << prod(*(result.L), *(result.U)) << std::endl;
//    char a =' ';
//    std::cin>>a;
    return 0;
}
