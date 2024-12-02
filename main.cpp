#include "set.hpp"

int main(){

    // Matrix a{
    //     {1}, 
    //     {1}, 
    //     {1},
    //     // {3}, 
    //     // {2}  

    // };
    // Matrix b{
    //     {3},
    //     {3}, 
    //     {3},
    //     // {2}, 
    //     // {2} 
    // };
    // Matrix c{
    //     {2},
    //     {2}, 
    //     {2},  
    //     // {3}, 
    //     // {3}
    // };
    
    // Matrix d{
    //     {0}, 
    //     {0}, 
    //     {1},
    //     // {3}, 
    //     // {3}

    // };
    // std::cout<<"DONE creating vectors\n";
    // std::cout<<"HERE IS AUGMENT\n"<<(a|b);
    // std::cout<<"DONE printing augment\n";


    // Matrix v1{
    //     {1}, 
    //     {2}, 
    //     {3},
    //     {4}
    // };
    // Matrix v2{
    //     {5}, 
    //     {6}, 
    //     {7}, 
    //     {8}
    // };
    Matrix f{
        {1}, 
        {2}, 
        {3}
    };
    Matrix g{
        {4}, 
        {5}, 
        {6}
    };

    // Matrix x{
    //     {1, 2, 1, 3, 4}, 
    //     {3, 6, 2, 6, 9}, 
    //     {-2, -4, 1, 1, -1}
    // };

    std::cout<<"BACK IN MAIN: f: "<<f;
    f|=g;
    std::cout<<"BACK IN MAIN: f: "<<f;

    // std::cout<<"______________________________________________________________\n";

    // std::cout<<"BACK IN MAIN: v1:"<<v1;
    // v1|=v2;
    // std::cout<<"BACK IN MAIN: v1:"<<v1;


    return 0;
}