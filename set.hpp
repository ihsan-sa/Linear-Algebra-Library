#ifndef SET_HPP
#define SET_HPP

#include "matrix.hpp"

class Node;
class Set;

class Node{
    Matrix value_;
    Node* p_next_;
public:
    Node(Matrix const &value, Node* p_next);
    Matrix value() const;
    Node* next() const;

    friend class Set;
};

class Set{
    Node* p_head_;
public:
//constructors and destructors
    Set(); 
    Set(std::initializer_list<Matrix> list, bool is_list);
    Set(Matrix const &m); //creates a set from the columns of a matrix
    Set(Set const &org); //copy construcyot
    Set(Set &&org); //move constructor

    Set &operator=(Set const &org); //copy assignment
    Set &operator=(Set &&org); //move assignment

    ~Set(); //destructor

//functions
    void clear(); //clears the set

    //inserting or removing
    int insert(Matrix const &m); //insert a matrix into a set
    int insert(std::initializer_list<Matrix> list); //insert a list of matrices into a set
    int insert(Set const &s); //inserts an entire set
    int remove(Matrix const &m); //removes a matrix from a set
    int remove(std::initializer_list<Matrix> list); //removes a list of matrices from a set

    Matrix convert_to_matrix() const; //returns the set but in matrix form

    //finding matrices in the set

    Node* find(Matrix const &m) const; //returns the address of the node if the matrix is found. Else, returns nullptr.
    bool contains(Matrix const &m) const; //will return true if the matrix m is contained in the set.
    bool in_span(Matrix const &m) const; //will return true if the matrix is in the span of s

    //Linear dependence and independence

    int dependencies() const; //returns the number of dependencies in a set
    int *dependencies_pos() const; //returns a DYNAMICALLY ALLOCATED ARRAY with the positions of the dependencies (cols in matrix)
    static int *dependencies_pos(Matrix const &m); //returns a DYNAMICALLY ALLOCATED ARRAY with the positions of the dependencies (cols in matrix)
    
    Set remove_dependencies() const; //removes all dependencies from the set and returns a new set
    
    int dim() const; //removes dependencies from a set S and then determines the dimension of SpanS
    
    //Systems
    static Set solve_system(Matrix const &A, Matrix const &b); //solves system and returns a set. First matrix in the set is the constant vector. Following matrices are the free variables.

    //fundamental subspaces associated with a matrix
    static Set col_space(Matrix const &m); //returns a basis for the column space of a matrix m
    static Set null_space(Matrix const &m); //returns a basis for the null space of a matrix m

    //eigenvalues and eigenvectors
    static Matrix evecs(Matrix const &m); //returns a matrix (P) with the evecs of a matrix.

    //printing a set
    friend std::ostream &operator<<(std::ostream &out, Set const &s);
};


#endif