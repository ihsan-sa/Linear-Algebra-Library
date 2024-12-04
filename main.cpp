#include "set.hpp"

int main(){

   Matrix A{
    {2, 1, 9}, 
    {0, 1, 2}, 
    {1, 0, 3}
   };
   Matrix b{
    {31}, 
    {8}, 
    {10}
   };

   std::cout<<Set::solve_system(A, b);

    return 0;
}