#include "set.hpp"

int main(){

//    Matrix A{
//     {2, 1, 9}, 
//     {0, 1, 2}, 
//     {1, 0, 3}
//    };
//    Matrix b{
//     {31}, 
//     {8}, 
//     {10}
//    };

    Matrix C{
        {1, 2, -3}, 
        {4, -5, 6},
        {-7, 8, 9}
    };

    Matrix D{
        {1,2,-1,3}, 
        {1,2,0,4}, 
        {0,0,0,3},
        {-1, 1, 2,1}
    };
    Matrix M{
        {1,2}, 
        {3,4}
    };

    Matrix A{
        {1, 1, 4, 6}, 
        {1, 2, 9, 12}, 
        {2, -1, 2, 3}
    };
    
    std::cout<<Set(Matrix{
        {-2, 0, 2}, 
        {1, 1, 0}, 
        {0, 2, 2}
    }).remove_dependencies();
   
    return 0;
}