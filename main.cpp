#include "set.hpp"

int main(){

    Matrix a{
        {1}, 
        {1}, 
        {1},
        // {3}, 
        // {2}  

    };
    Matrix b{
        {3},
        {3}, 
        {3},
        // {2}, 
        // {2} 
    };
    Matrix c{
        {2},
        {2}, 
        {2},  
        // {3}, 
        // {3}
    };
    
    Matrix d{
        {0}, 
        {0}, 
        {1},
        // {3}, 
        // {3}

    };
    // std::cout<<"DONE creating vectors\n";
    // std::cout<<"HERE IS AUGMENT\n"<<(a|b);
    // std::cout<<"DONE printing augment\n";



    Matrix e{
        {1}, 
        {0}
    };
    Matrix f{
        {1}, 
        {1}
    };
    Matrix g{
        {0}, 
        {1}
    };

    Matrix x{
        {1, 2, 1, 3, 4}, 
        {3, 6, 2, 6, 9}, 
        {-2, -4, 1, 1, -1}
    };


    Set null_x{Set::null_space(x)};
    Set col_x{Set::col_space(x)};
    std::cout<<null_x<<col_x;

    return 0;
}