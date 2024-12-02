#include "matrix.hpp"
#include <stdexcept>

#define DEBUG_STATUS false

//constructor and destructors

//question: Can you call an api func in a constructor

Matrix::Matrix() : matrix_{nullptr}, rows_{0}, cols_{0}{}
Matrix::Matrix(int rows, int cols) : matrix_{nullptr}, rows_{rows}, cols_{cols} {
    //dynamically allocate mem for the rows
    if(cols == -1){
        cols_ = rows_;
    }
    matrix_ = new float*[rows_];  //can i also call rows()?
    for(int m{0};  m <  rows; m++){
        matrix_[m] = new float[cols_]; // dynamically allocate mem for cols
        for(int n{0}; n < cols_; n++){
            matrix_[m][n] = 0; //set default value to 0
        }
    }
}   
Matrix::Matrix(std::initializer_list<std::initializer_list<float>> init_vals) : matrix_{nullptr}, rows_{0}, cols_{0} {
    rows_ = init_vals.size();
    for(std::initializer_list<float> row : init_vals){
        if(row.size() > cols_){
            cols_ = row.size();
        }
        //this way cols is equal to the biggest row
    }
    //initialize a new empty matrix
    matrix_ = new float*[rows_];
    for(int row{0}; row < rows_; row++){
        matrix_[row] = new float[cols_];
        for(int col{0}; col < cols_; col++){
            matrix_[row][col] = 0;
        }
    }

    //assign all the values
    int row_nbr = 0;
    for(std::initializer_list<float> row : init_vals){
        int col_nbr = 0;
        for(float col_val : row){
            matrix_[row_nbr][col_nbr] = col_val;
            col_nbr++;
        }
        row_nbr++;
    }
}
// Matrix::Matrix(std::initializer_list<float> init_vals) : matrix_{nullptr}, rows_{0}, cols_{1}{
//     rows_ = init_vals.size();
//     matrix_ = new float*[rows_];
//     int row{0};
//     for(float val : init_vals){
//         matrix_[row] = new float[cols_];
//         matrix_[row][0] = val;
//         row++;
//     }
// }

//copy constructor
Matrix::Matrix(Matrix const &org) : matrix_{nullptr}, rows_{org.rows()}, cols_{org.cols()}{
    //create a new matrix and set all the entries from the new matrix to the entries of the old one
    // std::cout<<"Calling copy constructor with"<<org;
    matrix_ = new float*[rows_];
    for(int row{0}; row < rows_; row++){
        matrix_[row] = new float[cols_];
        for(int col{0}; col < cols_; col++){
            matrix_[row][col] = org.matrix(row, col);
        }
    }
}
//copy assignment
Matrix &Matrix::operator=(Matrix const &org){
    // std::cout<<"Calling copy assignment with"<<org;

    clear();
    rows_ = org.rows();
    cols_ = org.cols();

    //create new matrix and assign everything
    matrix_ = new float*[rows()];
    for(int row{0}; row < rows(); row++){
        matrix_[row] = new float[cols()];
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] = org.matrix(row, col);
        }
    }
    return *this;
}
//move constructor
Matrix::Matrix(Matrix &&org) : matrix_{nullptr}, rows_{0}, cols_{0} {
    // std::cout<<"Calling move constructor with"<<org;

    //swap all member vars
    std::swap(matrix_, org.matrix_);
    std::swap(cols_, org.cols_);
    std::swap(rows_, org.rows_);
}
//move assignment
Matrix &Matrix::operator=(Matrix &&org){
    
    // std::cout<<"Calling move assignment with"<<org;

    //swap all member vars
    std::swap(matrix_, org.matrix_);
    std::swap(cols_, org.cols_);
    std::swap(rows_, org.rows_);

    return *this;
}
//destructor
Matrix::~Matrix(){
    // if(DEBUG_STATUS) std::cout<<"----------------------------------------------------------\n";
    if(DEBUG_STATUS) std::cout<<"\nCalling ~Matrix() on "<<*this;
    clear();
    if(DEBUG_STATUS) std::cout<<"Back in ~Matrix()\n\n";
}

//API
int Matrix::rows() const {
    return rows_;
}
int Matrix::cols() const {
    return cols_;
}
float Matrix::matrix(int m, int n) const {
   
    //safety: throw if trying to index out of range
    if(m >= rows() || n >= cols()){
        throw std::domain_error{ //SHOULD THIS BE DOMAIN ERR
            "Attempted to access index (" + std::to_string(m) + "," + std::to_string(n) + ") in matrix with dimensions " + std::to_string(rows()) + " x " + std::to_string(cols()) +". Index out of range!"
        };
        return -1; //is there something better to return?
    }
    return matrix_[m][n];
}

//Operators

std::ostream &operator<<(std::ostream &out, Matrix const &m){

    out<<std::endl;
    for(int row{0}; row < m.rows(); row++){
        for(int col{0}; col < m.cols(); col++){
            out<<m.matrix(row, col)<<" ";
        }
        out<<std::endl;
    }
    out<<std::endl;
    return out;
}

Matrix operator*(float const &scalar, Matrix const &m){
    //multiply each entry by the scalar
    Matrix tmp{m};
    for(int row{0}; row < tmp.rows();  row++){
        for(int col{0}; col < tmp.cols(); col++){
            tmp.matrix_[row][col] *= scalar;
        }
    }
    return tmp;
}
Matrix &Matrix::operator+=(Matrix const &rhs){
    if(!is_same_dim(rhs)){
        throw std::domain_error{
            "Error, dimensions do not agree for addition."
        };
    }
    //for each entry, increment by the corresponding entry in rhs
    for(int row{0}; row < rows();  row++){
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] += rhs.matrix(row, col);
        }
    }
    return *this;
}
Matrix &Matrix::operator-=(Matrix const &rhs){
    if(!is_same_dim(rhs)){
        throw std::domain_error{
            "Error, dimensions do not agree for subtraction."
        };
    }
    //for each entry, increment by the corresponding entry in rhs
    for(int row{0}; row < rows();  row++){
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] -= rhs.matrix(row, col);
        }
    }
    return *this;
}
Matrix &Matrix::operator*=(Matrix const &rhs){
    if(!is_mult_allowed(rhs)){
        throw std::domain_error{
            "Dimensions do not agree for matrix multiplication."
        };
    }
    //assign a tmp matrix
    Matrix tmp = *this;
    //clear *this
    clear();
    //redefine the dimensions. nbr of rows stays the same
    cols_ = rhs.cols(); 
    rows_ = tmp.rows();

    matrix_ = new float*[rows()];
    for(int row{0}; row <  rows(); row++){
        matrix_[row] = new float[cols()];
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] = tmp.vrow(row)%rhs[col];
        }
    }

    return *this;
}
Matrix &Matrix::operator*=(float const &scalar){
    //multiply each entry by the scalar
    for(int row{0}; row < rows();  row++){
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] *= scalar;
        }
    }
    return *this;
}
Matrix &Matrix::operator/=(Matrix const &rhs){
    Matrix rref_aug{(*this|rhs).rref()};
    int cols_in_org{cols()};

    clear();
    *this = rhs;

    for(int row{0}; row < rows(); row++){
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] = rref_aug.matrix(row, col + cols_in_org);
        }
    }

    return *this;
}
Matrix &Matrix::operator|=(Matrix const &rhs){
    if(DEBUG_STATUS) std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    if(DEBUG_STATUS) std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    if(DEBUG_STATUS) std::cout<<"\nIn |=\n";
    if(DEBUG_STATUS) std::cout<<*this<<rhs<<"\t*this dim: "<<rows()<<"x"<<cols()<<" rhs: "<<rhs.rows()<<"x"<<rhs.cols()<<"\n";
    if(rows() != rhs.rows()){
        throw std::domain_error{
            "Number of rows does not agree for augment."
        };
    }
    //create temp matrix
    Matrix tmp{*this};
    if(DEBUG_STATUS) std::cout<<"\t*this.matrix_: "<<matrix_<<" tmp.matrix_: "<<tmp.matrix_<<"\n";
    //clear current matrix
    if(DEBUG_STATUS) std::cout<<"\tAbout to clear *this\n";
    clear();
    if(DEBUG_STATUS) std::cout<<"\tCleared *this\n";
    //reassign matrix
    if(DEBUG_STATUS) std::cout<<"\t*this.matrix_: "<<matrix_<<" tmp.matrix_: "<<tmp.matrix_<<"\n";
    cols_ = tmp.cols() + rhs.cols();
    rows_ = tmp.rows();
    matrix_ = new float*[cols()];
    if(DEBUG_STATUS) std::cout<<"\tReassigned matrix_: "<<matrix_<<"\n";


    //now reassign temp
    for(int row{0}; row < rows(); row++){
        matrix_[row] = new float[cols()];
        for(int col{0}; col < tmp.cols(); col++){
            matrix_[row][col] = tmp.matrix(row, col);
        }
    }
    if(DEBUG_STATUS) std::cout<<*this;
    //now perform augment
    for(int row{0}; row < rows(); row++){
        for(int col{0}; col < rhs.cols(); col++){
            matrix_[row][tmp.cols() + col] = rhs.matrix(row, col);
        }
    }
    if(DEBUG_STATUS) std::cout<<*this;
    if(DEBUG_STATUS) std::cout<<"\t*this.matrix_: "<<matrix_<<" tmp.matrix_: "<<tmp.matrix_<<" *matrix_: "<<*tmp.matrix_<<"\n";
    if(DEBUG_STATUS) std::cout<<"DONE |=\n";
    if(DEBUG_STATUS) std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    if(DEBUG_STATUS) std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";


    return *this;
}
float Matrix::operator%(Matrix const &rhs) const{
    if(!is_vector() || !rhs.is_vector()){
        throw std::domain_error{
            "Cannot perform dot product on non-vector matrices."
        };
    }
    if(!is_same_dim(rhs)){
        throw std::domain_error{
            "Vectors are not in the same space." //is Rn a space?
        };
    }
    //perform dot product
    float dot_product = 0;
    for(int row{0}; row < rows(); row++){
        dot_product += matrix(row, 0) * rhs.matrix(row, 0);
    }
    return dot_product;
    
}
Matrix Matrix::operator[](int col) const{
    Matrix tmp(rows(), 1); //create a new vector
    for(int row{0}; row < rows(); row++){
        tmp.matrix_[row][0] = matrix(row, col);
    }
    return tmp;
}

Matrix Matrix::operator+(Matrix const &rhs) const{
    Matrix tmp{*this};
    tmp += rhs;
    return tmp;
}
Matrix Matrix::operator-(Matrix const &rhs) const{
    Matrix tmp{*this};
    tmp -= rhs;
    return tmp;
}
Matrix Matrix::operator*(Matrix const &rhs) const{
    Matrix tmp{*this};
    tmp *= rhs;
    return tmp;
}
Matrix Matrix::operator/(Matrix const &rhs) const{
    Matrix tmp{*this};
    tmp /= rhs;
    return tmp;
}
Matrix Matrix::operator|(Matrix const &rhs) const{
    Matrix tmp{*this};
    tmp |= rhs;
    return tmp;
}
Matrix Matrix::operator*(float const &scalar) const{
    Matrix tmp{*this};
    tmp *= scalar;
    return tmp;
}
//Bool ops

bool Matrix::is_vector() const {
    return (cols() == 1) ? true : false;
}
bool Matrix::is_square() const { 
    return (rows() == cols()) ? true : false;
}
bool Matrix::operator==(Matrix const &m) const {
    if(!is_same_dim(m)) return false; //if dimensions are not the same, return 0;

    for(int row{0}; row < rows(); row++){
        for(int col{0}; col < cols(); col++){
            if(matrix(row, col) != m.matrix(row, col)) return false; //if an entry doesn't match, return 0;
        }
    }
    //if we get here, both matrices are the same
    return true;
}
bool Matrix::operator!=(Matrix const &m) const {
    return !(*this == m); //return the negation of the comparison
}
bool Matrix::is_mult_allowed(Matrix const &m) const {
    return (cols() == m.rows()) ? true : false;
}
bool Matrix::is_same_dim(Matrix const &m) const {
    return ((rows() == m.rows()) && (cols() == m.cols())) ? true : false;
}
bool Matrix::is_vector(Matrix const &m){
    return m.is_vector();
}
bool Matrix::is_square(Matrix const &m){
    return m.is_square();
}
bool Matrix::is_mult_allowed(Matrix const &m1, Matrix const &m2){
    return m1.is_mult_allowed(m2);
}
bool Matrix::is_same_dim(Matrix const &m1, Matrix const &m2){
    return m1.is_same_dim(m2);
}

//operations
void Matrix::clear(){
    if(DEBUG_STATUS) std::cout<<"============================================================================================\n";
    if(DEBUG_STATUS) std::cout<<"=  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n";

    if(DEBUG_STATUS) std::cout<<"In clear\n";
    if(matrix_ == nullptr){
        if(DEBUG_STATUS) std::cout<<"\tMatrix already cleared.";
        return;
    }
    if(DEBUG_STATUS) std::cout<<"\tGoing to delete rows\n\n";
    for(int row{0}; row < rows(); row++){
        for(int col{0}; col < cols(); col++){
            matrix_[row][col] = 0;
        }
        if(DEBUG_STATUS) std::cout<<"\t\tabout to delete row: "<<row<<" "<<matrix_[row]<<"\n";
        delete[] matrix_[row];
        if(DEBUG_STATUS) std::cout<<"\t\tdeleted row: "<<row<<"\n";
        matrix_[row] = nullptr;
    }
    if(DEBUG_STATUS) std::cout<<"\n\tabout to delete matrix: "<<matrix_<<" "<<*matrix_<<"\n";
    delete[] matrix_;
    if(DEBUG_STATUS) std::cout<<"\tdeleted matrix\n";
    if(DEBUG_STATUS) std::cout<<"Done clear\n";
    if(DEBUG_STATUS) std::cout<<"=  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n";  
    if(DEBUG_STATUS) std::cout<<"============================================================================================\n";


    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
}
Matrix Matrix::row(int row) const{
    //safety
    if(row >= rows()){
        throw std::domain_error{
            "Row out of range of matrix."
        };
    }
    
    Matrix tmp{1, cols()};
    for(int col{0}; col < cols(); col++){
        tmp.matrix_[0][col] = matrix(row, col);
    }
    return tmp;
}
Matrix Matrix::vrow(int row) const{
    //safety
    if(row >= rows()){
        throw std::domain_error{
            "Row out of range of matrix."
        };
    }
    
    Matrix tmp{cols(), 1};
    for(int col{0}; col < cols(); col++){
        tmp.matrix_[col][0] = matrix(row, col);
    }
    return tmp;
}
Matrix Matrix::eye(int rows, int cols){
    if(cols == -1) cols = rows; //is there a way of making this safer? The comparison could accidentally become an assignment and cause error.
    Matrix tmp{rows, cols};
    
    for(int row{0}, col{0}; row < tmp.rows() && col < tmp.cols(); row++, col++){
        tmp.matrix_[row][col] = 1;
    }
    return tmp;
}
Matrix Matrix::ref() const{

    float constexpr TOLERANCE = 10e-3;

    Matrix ref{*this};
    if(ref.rows() == 0) return ref;

    //first, we start with column 0

    for(int col{0}, row{0}; col < ref.cols() && row < ref.rows(); col++, row++){
        // std::cout<<"Current iteration: ("<<row<<","<<col<<")\n"<<ref;
        int iter_row{row};
        bool no_nonzero_found{false};
        
        while(std::abs(ref.matrix(iter_row, col)) < TOLERANCE){
            // std::cout<<"\tLooking for nonzero entry... on row: "<<iter_row<<"\n";
            iter_row++;
            if(iter_row == ref.rows()){
                // std::cout<<"\tNo nonzero entry in this column... going to next.\n";
                row--; //if we get to the end of the column and there was no nonzero entry, we move to the next column
                no_nonzero_found = true;
                break;
            }
        }
        if(no_nonzero_found) continue;

        //yay, we found a nonzero entry! Now we swap
        // std::cout<<"\tSwapping\n";
        ref.swap_rows(row, iter_row);

        //now, for every row below row, we subtract R1 according to algorithm

        float row_le = ref.matrix(row, col);
        // std::cout<<"\tLeading entry: "<<row_le<<"\n";

        for(int sub_row{row + 1}; sub_row < ref.rows(); sub_row++){

            // std::cout<<"\t\tIterating thru rows... on row: "<<sub_row<<"\n";
            float sub_row_le = ref.matrix(sub_row, col); //get the leading entry
            float factor = sub_row_le/row_le;

            for(int sub_col{col}; sub_col < ref.cols(); sub_col++){
                // std::cout<<"\t\t\tIterating thru cols... on col: "<<sub_col<<"\n";
                ref.matrix_[sub_row][sub_col] -= ref.matrix(row, sub_col) * factor;
                //if the entry is smaller than the tolerance, we set it to zero
                if(std::abs(ref.matrix(sub_row, sub_col)) < TOLERANCE){
                    ref.matrix_[sub_row][sub_col] = 0.0f;
                }

            }

            // std::cout<<"\t\tDone sub row: "<<sub_row<<"\n";

        }

    }

    // std::cout<<"Done ref.\n"<<ref;
    return ref;

    
}
Matrix Matrix::rref() const{
    Matrix rref = ref();
    // std::cout<<"In rref()\n";
    // std::cout<<rref;
    constexpr int TOLERANCE = 10e-3;

    //first we set all of the leading entries to 1
    // std::cout<<"\tAbout to normalize all leading entries.\n\n";

    for(int row{0}; row < rref.rows(); row++){
        // std::cout<<"\n\t\tNormalizing... row: "<<row<<"\n";

        bool leading_entry_found = true;
        int col{0};
        for(; rref.matrix(row, col) == 0; col++){
            // std::cout<<"\t\t\tLooking for leading entry... col: "<<col<<"\n";

            if(col == rref.cols() - 1){
                //we are on the last entry so break and go to the next row
                leading_entry_found = false;
                // std::cout<<"\t\tLeading entry not found. Continuing\n";
                break;
            }
        }

        if(!leading_entry_found) continue;

        float leading_entry = rref.matrix(row, col);
        // std::cout<<"\t\tLeading entry found: "<<leading_entry<<"\n";

        //now we loop through all the elements in that row and divide them by the leading entry
        // std::cout<<"\t\tAbout to set all entries...\n";
        for(int col_loop{0}; col_loop < rref.cols(); col_loop++){
            rref.matrix_[row][col_loop] /= leading_entry;

            if(std::abs(rref.matrix(row, col_loop)) < TOLERANCE) rref.matrix_[row][col_loop] = 0.0f; //account for deviations
            if(rref.matrix(row, col_loop) == -0) rref.matrix_[row][col_loop] = 0.0f; //acount for -0

            // std::cout<<"\t\t\tSet col: "<<col_loop<<" to: "<<rref.matrix(row, col_loop)<<"\n";

        }

    }


    // std::cout<<"\n\tAll leading entries have been normalized!\n";

    //Now we need to start at the bottom, find the leading entry and then subtract the above rows

    for(int row{rref.rows() - 1}; row >= 0; row--){
        // std::cout<<"\n\t\tPerforming rref.. row: "<<row<<"\n";

        bool leading_entry_found = true;
        int col{0};
        for(; rref.matrix(row, col) == 0; col++){
            // std::cout<<"\t\t\tLooking for leading entry... col: "<<col<<"\n";

            if(col == rref.cols() - 1){
                //we are on the last entry so break and go to the next row
                leading_entry_found = false;
                // std::cout<<"\t\tLeading entry not found. Continuing\n";
                break;
            }
        }

        if(!leading_entry_found) continue;

        float leading_entry = rref.matrix(row, col);
        // std::cout<<"\t\tLeading entry found: "<<leading_entry<<" at col: "<<col<<"\n";

        //now for each row above it, we sub

        for(int row_above{row - 1}; row_above >= 0; row_above--){
            float entry_above_le = rref.matrix(row_above, col);
            float mult_factor = entry_above_le/leading_entry;

            for(int col_loop{0}; col_loop < rref.cols(); col_loop++){

                rref.matrix_[row_above][col_loop] -= rref.matrix(row, col_loop) * mult_factor;

                if(std::abs(rref.matrix(row_above, col_loop)) < TOLERANCE) rref.matrix_[row_above][col_loop] = 0.0f; //account for deviations
                if(rref.matrix(row_above, col_loop) == -0) rref.matrix_[row_above][col_loop] = 0.0f; //acount for -0

                // std::cout<<"\t\t\tSet col: "<<col_loop<<" to: "<<rref.matrix(row_above, col_loop)<<"\n";

            }
        }
    }

    // std::cout<<"Done rref!"<<rref;

    return rref;

}
Matrix &Matrix::swap_rows(int r1, int r2){
    if((r1 > rows() || r2 > rows()) || (r1 < 0 || r2 < 0)){
        throw std::domain_error{
            "Selected rows for row swap are not within range."
        };
    }
    Matrix row_a = row(r1);
    Matrix row_b = row(r2);

    for(int col{0}; col < cols(); col++){
        //swap rows
        matrix_[r1][col] = row_b.matrix(0, col);
        matrix_[r2][col] = row_a.matrix(0, col);
    }

    return *this;

}
int Matrix::rank() const{

    Matrix ref_matrix{ref()};

    int rank{0};
    for(int row{0}; row < rows(); row++){
        for(int col{0}; col < cols(); col++){
            if(ref_matrix.matrix(row, col) != 0){
                rank++;
                break;
            }
        }
    }
    return rank;
}
int Matrix::nullity() const{
    return cols() - rank();
}
Matrix Matrix::inv() const{
    if(!is_square()){
        throw std::domain_error{
            "Non square matrix cannot be inverted."
        };
    }
    return Matrix{(*this)/(Matrix::eye(cols()))};
}
Matrix Matrix::transpose() const{
    Matrix tmp{cols(), rows()};

    for(int row{0}; row < tmp.rows(); row++){
        for(int col{0}; col < tmp.cols(); col++){
            tmp.matrix_[row][col] = matrix(col, row);
        }
    }
    return tmp;
}
Matrix Matrix::e(int n, int i){
    if(i > n){
        throw std::domain_error{
            "index out of range."
        };
    }
    Matrix tmp{n, 1}; //create vector
    tmp.matrix_[i - 1][0] = 1;
    return tmp;
}
float Matrix::get_le_val(int row) const {
    if(row >= rows()){
        throw std::domain_error{
            "Index out of range - leading entry could not be found."
        };
    }

    bool leading_entry_found = true;
    int col{0};

    for(; matrix(row, col) == 0; col++){
        std::cout<<"\t\t\tLooking for leading entry... col: "<<col<<"\n";

        if(col == cols() - 1){
            //we are on the last entry so break and go to the next row
            leading_entry_found = false;
            std::cout<<"\t\tLeading entry not found.\n";
            break;
        }
    }

    if(!leading_entry_found) return 0;

    float leading_entry = matrix(row, col);
    std::cout<<"\t\tLeading entry found: "<<leading_entry<<" at col: "<<col<<"\n";

    return leading_entry;
}
int Matrix::get_le_col(int row) const{
    if(row >= rows()){
        throw std::domain_error{
            "Index out of range - leading entry could not be found."
        };
    }

    bool leading_entry_found = true;
    int col{0};

    for(; matrix(row, col) == 0; col++){

        if(col == cols() - 1){
            //we are on the last entry so break and go to the next row
            leading_entry_found = false;
            break;
        }
    }

    if(!leading_entry_found) return -1;


    return col;
}




























//old ref
// Matrix Matrix::ref() const{
//     Matrix ref{*this};
//     constexpr float TOLERANCE = 10e-3;
//     int min_m_n = std::min(ref.rows(), ref.cols()); //minimum btw the number of rows and the number of columns

//     //now we loop through each row until the min_m_n and find the leading entry

//     for(int row{0}, offset{0}; row < min_m_n; row++){
//         float leading_entry = ref.matrix(row, row + offset);

//         //if the leading entry is zero, we loop through the rest of the rows to find if there is a row with nonzero leading entry.
//         //if we find one, we swap them
//         if(leading_entry < TOLERANCE){

//             for(int new_row{row + 1}; new_row < ref.rows(); new_row++){

//                 if(ref.matrix(new_row, row + offset) != 0){

//                     leading_entry = ref.matrix(new_row, row + offset);
//                     ref.swap_rows(row, new_row);
//                     break;
//                 }
//             }

//         }
        
//         //is this necessary
//         if(std::abs(leading_entry) < TOLERANCE){
//             // offset++;
//             // min_m_n--;
//             // row--;
//             continue;
    
//         }


//         //now, for each row below it, we subtract the first row divided by leading entry times the first entry of the row

//         for(int next_row{row + 1}; next_row < ref.rows(); next_row++){
//             float row_first_entry = ref.matrix(next_row, row + offset);
//             float factor = row_first_entry/leading_entry;
//             //now, for each entry in that row, we subtract the corresponding entry in the first row as explained by the algorithm

//             for(int col{row}; col < ref.cols(); col++){
//                 ref.matrix_[next_row][col] -= ref.matrix(row, col) * factor;
//                 if(std::abs(ref.matrix(next_row, col)) < TOLERANCE){
//                     ref.matrix_[next_row][col] = 0.0f;
//                 }

//             }

//         }
//     }

//     return ref;
// }


//OLD RREF

// Matrix Matrix::rref() const{

//     Matrix rref{ref()};
//     //start at bottom and work our way up
//     for(int row{std::min(cols(), rows())-1}, entry{0};  row > 0; row--, entry++){
        
//         float last_entry = rref.row(row).matrix_[0][std::min(cols(), rows()) - 1 - entry]; 
 
//         for(int k{row - 1}; k >= 0; k--){
//             float row_le = rref.row(k).matrix_[0][std::min(cols(), rows())- 1 - entry];
//             // std::cout<<" next row: "<<k<<" row_fe: "<<row_le<<std::endl;
//             // subtract prev row * first entry next row from next row

//             for(int i{0}; i < cols(); i++){
//                 // std::cout<<"   "<<matrix_rref.matrix_[k][i]<<" - "<< matrix_rref.matrix_[row][i] * ((last_entry == 0) ? 0 : (row_le/last_entry))<< " = "<<matrix_rref.matrix_[k][i] - matrix_rref.matrix_[row][i] *  ((last_entry == 0) ? 0 : (row_le/last_entry))<<std::endl;
//                 rref.matrix_[k][i] -= rref.matrix(row, i) * ((last_entry == 0) ? 0 : (row_le/last_entry));
                
//             }
//         }
        

//         // std::cout<<"\n";
//         // matrix_rref.print();
//         // std::cout<<std::endl;
//     }

//     for(int i{0}; i < std::min(cols(), rows()); i++){
//         float leading_entry = 1;
//         for(int j{0}; j<cols(); j++){
//             if(rref.matrix_[i][j] != 0){
//                 leading_entry = rref.matrix(i, j);
//                 break;
//             }
//         }
//         for(int j{0}; j < rref.cols();j++){
//             // std::cout<<matrix_rref.matrix_[i][j] <<"/"<< leading_entry<<" = "<<matrix_rref.matrix_[i][j] / leading_entry<<std::endl;
//             rref.matrix_[i][j] /= leading_entry;
            
//             if(rref.matrix_[i][j] == -0){
//                 rref.matrix_[i][j] = 0;
//             }

//         }
        

//         // std::cout<<std::endl;
//     }
//     return rref;   
// }

