#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <initializer_list>
// #include "set.hpp"

class Matrix;
class Set;

class Matrix{
    float **matrix_;
    int rows_;
    int cols_;
 
public:

//constructors and destructor

    //construtors
    Matrix(); //default constructor -- DONE
    Matrix(int rows, int cols = -1); //init an mxn matrix full of zeros. If cols = -1, then creates a square matrix -- DONE
    Matrix(std::initializer_list<std::initializer_list<float>> init_vals); //initialize a matrix using a std initializer list. -- DONE
    // Matrix(std::initializer_list<float> init_vals); //initialize a vector

    //destructor
    ~Matrix(); // -- DONE

    //Copy constructor and copy assignment
    Matrix(Matrix const &org); // -- DONE
    Matrix &operator=(Matrix const &org); // -- DONE

    //Move constructor and move assignment
    Matrix(Matrix &&org); // -- DONE
    Matrix &operator=(Matrix &&org); // -- DONE

//API

    int rows() const; //returns the number of rows -- DONE
    int cols() const; //return the number of cols  -- DONE
    float matrix(int row, int col) const; //returns the entry at (row, col) -- DONE

//question: Keep in mind, is there any function which should return a constant matrix or a constant reference
//probaby not because then we can't modify in the future...
//Basic operations

    //increments
    Matrix &operator+=(Matrix const &rhs); //add rhs to *this and return *this by reference. *this is modified -- DONE
    Matrix &operator-=(Matrix const &rhs); //subtract rhs from *this and return *this by reference. *this is modified -- DONE
    Matrix &operator*=(Matrix const &rhs); //matrix multiplication between *this and rhs and return *this by reference. *this is modified and dimensions change -- DONE
    Matrix &operator*=(float const &scalar); //scalar mult. modifies and returns *this. &scalar is not necessary but could help save mem in long run -- DONE
    Matrix &operator/=(Matrix const &rhs); //solve the system *this * x = rhs and reassign result x to *this. Return *this by reference. *this is modified.
    Matrix &operator|=(Matrix const &rhs); //augment rhs onto *this and return *this by ref. *this is modified. -- DONE

    //operations
    Matrix operator+(Matrix const &rhs) const; //addition: returns a brand new matrix and does not change the matrix it is called on -- DONE
    Matrix operator-(Matrix const &rhs) const; //subtraction: returns a brand new matrix and does not change the matrix it is called on -- DONE
    Matrix operator*(Matrix const &rhs) const; //matrix multiplication: returns a brand new matrix and does not change the matrix it is called on -- DONE
    Matrix operator*(float const &scalar) const; // scalar multiplicaion: returns brand new matrix and does not change *thos.
    Matrix operator/(Matrix const &rhs) const; //solve the system *this * x = rhs and return brand new matrix x. *this is not modified.
    Matrix operator|(Matrix const &rhs) const; //returns the brand new augmented matrix *this|rhs. *this is not modified -- DONE
    Matrix operator[](int col) const; //returns the column at index col in a matrix as a brand new matrix (vector). *this is not modified. -- DONE
    float operator%(Matrix const &rhs) const; //returns the dot product between two vectors. *this is not modified. -- DONE

    //question: is this a good way of using friend?
    friend std::ostream &operator<<(std::ostream &out, Matrix const &m); //cout  -- DONE
    friend Matrix operator*(float const &scalar, Matrix const &m); //scalar mult with scalar in front. Creates a brand new matrix -- DONE

//Bool checks

    //operations on a matrix
    bool is_vector() const; //returns true if *this is a vector. Does not modify *this. -- DONE
    bool is_square() const; //returns true if *this is a square matrix. -- DONE
    bool operator==(Matrix const &m) const; //returns true if both matrices are equal -- DONE
    bool operator!=(Matrix const &m) const; //returns true if both matrices not are equal -- DONE
    bool is_mult_allowed(Matrix const &m) const; //returns true if *this * m is allowed -- DONE
    bool is_same_dim(Matrix const &m) const; //returns true if the dimensions of *this and m are the same -- DONE
    bool is_row_zero(int row) const; //returns true if the row is all zeros

    //static functions for the same operations
    static bool is_vector(Matrix const &m); //returns true if *this is a vector. Does not modify *this. -- DONE
    static bool is_square(Matrix const &m); //returns true if *this is a square matrix. -- DONE
    static bool is_mult_allowed(Matrix const &m1, Matrix const &m2); //returns true if *this * m is allowed -- DONE
    static bool is_same_dim(Matrix const &m1, Matrix const &m2); //returns true if the dimensions of *this and m are the same -- DONE
    

//Matrix functions

    //these functions return a brand new matrix. Perhaps it could eventually be preferable to return *this by reference and bring that matrix to ref
    Matrix ref() const; //returns a brand new matrix which is the ref of *this. Does not modify *this.
    Matrix rref() const; //returns a brand new matrix which is the rref of *this. Does not modify *this.
    Matrix inv() const; //returns a brand new matrix which is the inverse of *this. Does not modify *this.
    Matrix transpose() const; //returns a brand new matrix which is the transpose of *this. Does not modify *this.
    float det() const; //returns the determinant of *this. Does not modify *this.
    int rank() const; //returns the rank of *this. Does not modify *this.
    int nullity() const; //returns the nullity of *this. Does not modify *this.
    Matrix row(int row) const; //returns a brand new matrix corresponding to the row at index row. -- DONE
    Matrix vrow(int row) const; //returns a brand new vector corresponding to the transpose of the row at index row. -- DONE
    void clear(); //sets all entries to zero and deletes everything -- DONE
    Matrix &swap_rows(int r1, int r2); //swaps the rows and returns the matrix by reference
    float get_le_val(int row) const; //returns the leading entry at row row. Returns 0 if leading entry is not found
    int get_le_col(int row) const; //returns the column of the leading entry
    Matrix remove_zero_rows() const; //returns a new matrix which is a copy of *this with no zero rows
    Matrix remove_col(int col) const; //returns a new matrix which is a copy of *this with col col removed.
    Matrix remove_row(int row) const; //returns a new matrix which is a copy of *this with row row removed.

    //static functions

    static Matrix eye(int rows, int cols = -1); //retuens the rowsxcols identity matrix. If cols is not specified, then the matrix is square.  -- DONE
    static Matrix e(int n, int i); //returns the unit vector in Rn with the ith entry as one. Careful, does not use zero-base indexing

    friend class Set;
};


#endif