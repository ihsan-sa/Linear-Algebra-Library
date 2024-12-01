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

    Set &operator=(Set const &org);
    Set &operator=(Set &&org);

    ~Set();
//functions
    void clear();

    int insert(Matrix const &m);
    int insert(std::initializer_list<Matrix> list);
    int remove(Matrix const &m);
    int remove(std::initializer_list<Matrix> list);
    Matrix convert_to_matrix() const; //returns the set but in matrix form

    Node* find(Matrix const &m) const;

    bool contains(Matrix const &m) const;
    bool in_span(Matrix const &m) const;
    int dependencies() const; //returns the number of dependencies
    int *dependencies_pos() const; //returns a DYNAMICALLY ALLOCATED ARRAY with the positions of the dependencies
    Set remove_dependencies() const; //removes all dependencies from the set and returns a new set
    int dim() const; //removes dependencies from a set S and then determines the dimension of SpanS

    //fundamental subspaces associated with a matrix
    static Set col_space(Matrix const &m); //returns a basis for the column space of a matrix m
    static Set null_space(Matrix const &m); //returns a basis for the null space of a matrix m

    friend std::ostream &operator<<(std::ostream &out, Set const &s);
};




#endif