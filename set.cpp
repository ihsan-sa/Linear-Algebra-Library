#include "set.hpp"

//Node
Node::Node(Matrix const &value, Node* p_next) : value_{value}, p_next_{p_next}{

}

Node* Node::next() const {
    return p_next_;
}
Matrix Node::value() const {
    return value_;
}

//Set
Set::Set() : p_head_{nullptr}{
    
}
Set::Set(std::initializer_list<Matrix> list, bool is_list) : p_head_{nullptr}{
    for(Matrix m : list){
        insert(m);
    }
}
Set::Set(Set const &org) : p_head_{nullptr}{

    //loop through org
    for(Node* p_ptr{org.p_head_}; p_ptr != nullptr; p_ptr = p_ptr->next()){
        insert(p_ptr->value());
    }

}   
Set::Set(Set &&org) : p_head_{nullptr}{
    std::swap(p_head_, org.p_head_);
}
Set::Set(Matrix const &m) : p_head_{nullptr}{
    for(int col{0}; col < m.cols(); col++){
        insert(m[col]);
    }
}

Set &Set::operator=(Set const &org){
    if(this == &org) return *this;
    clear();
    for(Node* p_ptr{org.p_head_}; p_ptr != nullptr; p_ptr = p_ptr->next()){
        insert(p_ptr->value());
    }
    return *this;
}
Set &Set::operator=(Set &&org){
    std::swap(p_head_, org.p_head_);
    return *this;
}
// Move assignment


Set::~Set(){
    // std::cout<<"Calling ~Set()\n";
    clear();
}

void Set::clear(){
    for(Node* p_ptr{p_head_}; p_ptr != nullptr;){

        Node* p_next = p_ptr->next();
        delete p_ptr;
        p_ptr = p_next;
    }

    p_head_ = nullptr;
}
Node* Set::find(Matrix const &m) const{
    for(Node* p_ptr{p_head_}; p_ptr != nullptr; p_ptr = p_ptr->next()){
        if(p_ptr->value() == m){
            return p_ptr;
        }
    }
    return nullptr;
} 

int Set::insert(Matrix const &m){
    if(find(m) != nullptr) return 0;

    //we want to insert at the back
    Node* p_prev{nullptr};
    for(Node* p_ptr{p_head_}; p_ptr != nullptr; p_prev = p_ptr, p_ptr = p_ptr->p_next_); //loop until we get to the end

    if(p_prev == nullptr) p_head_ = new Node{m, nullptr};
    else p_prev->p_next_ = new Node{m, nullptr};

    return 1;
}
int Set::insert(std::initializer_list<Matrix> list){
    int inserted{0};
    for(Matrix m : list){
        inserted += insert(m);
    }
    return inserted;
}
int Set::remove(Matrix const &m){
    Node* p_prev{nullptr};
    for(Node* p_ptr{p_head_}; p_ptr != nullptr; p_prev = p_ptr, p_ptr = p_ptr->next()){
        if(p_ptr == find(m)){
            if(p_prev == nullptr){
                p_head_ = p_ptr->next();
                delete p_ptr;
                return 1;
            }
            p_prev->p_next_ = p_ptr->next();
            delete p_ptr;
            return 1;
        }
    }
    return 1;
}
int Set::remove(std::initializer_list<Matrix> list){
    int removed{0};
    for(Matrix m : list){
        removed += remove(m);
    }
    return removed;
}


bool Set::contains(Matrix const &m) const{
    return (find(m) == nullptr) ? false : true;
}
Matrix Set::convert_to_matrix() const{
//first create a matrix with all the vectors in the set
    if(p_head_ == nullptr) {
        throw std::domain_error{
            "Set empty."
        };
    } //update the dimensions

    Matrix spanning_set{p_head_->value()};

    // std::cout<<"Updated rows\n\n";

    for(Node* p_ptr{p_head_->next()}; p_ptr != nullptr; p_ptr = p_ptr->next()){
        
        // std::cout<<"Looping through set: "<<p_ptr<<std::endl;

        if(!p_ptr->value().is_vector()){
            throw std::domain_error{
                "Cannot create matrix with non vector quantities. In in_span()"
            };
        }
        if(p_ptr->value().rows() != spanning_set.rows()){
            throw std::domain_error{
                "Dimensions of vectors do not all agree."
            };
        }
        // std::cout<<"Augmenting on: "<<p_ptr->value();
        spanning_set |= p_ptr->value(); // augment on the matrix
        // std::cout<<"Done augment\n\n";
    }
    // std::cout<<"All done. "<<spanning_set;
    // std::cout<<spanning_set;
    return spanning_set;
}
bool Set::in_span(Matrix const &m) const{
    Matrix spanning_set{convert_to_matrix()};

    if(!m.is_vector() || m.rows() != spanning_set.rows()){
        throw std::domain_error{
            "Dimensions of m do not agree with set."
        };
    }
    //now we augment on the matrix
    Matrix augmented{spanning_set | m};
    // now we compare rank

    return (spanning_set.rank() == augmented.rank()) ? true : false;
}
int Set::dependencies() const{
    //in order to determine whether or not a system is LI, we need to solve the homogeneous system. 
    //if there is only one solution, then the system is LI.

    Matrix set{convert_to_matrix()};
    return set.cols() - set.rank();
}
std::ostream &operator<<(std::ostream &out, Set const &s){
    out<<"Set: \n";
    for(Node* p_ptr{s.p_head_}; p_ptr != nullptr; p_ptr = p_ptr->next()){
        out<<p_ptr->value();
    }
    return out;
}
Set Set::remove_dependencies() const{
    // std::cout<<"In remove_dependencies()\n";
    Set li_set{*this};
    Matrix ld_matrix{li_set.convert_to_matrix()};
    Matrix li_matrix{li_set.convert_to_matrix().rref()};
    // std::cout<<"\tli_matrix: ";
    // std::cout<<li_matrix;

    int row_offset{0};

    for(int col{0}; col < li_matrix.cols(); col++){
        // std::cout<<"\tLooking through matrix... row: "<<col+row_offset<<" col: "<<col<<" col - offset: "<<col - row_offset<<" le: "<<li_matrix.matrix(col + row_offset, col)<<"\n"; 

        if(col + row_offset == li_matrix.rows() - 1){
            // std::cout<<"Removing remaining entries\n";
            for(int remaining_col{col + 1}; remaining_col < li_matrix.cols(); remaining_col++){
                li_set.remove(ld_matrix[remaining_col]); //we remove all the next cols
            }
            break;
        }

        if(li_matrix.matrix(col + row_offset, col) != 1){
            // std::cout<<"\tAbout to remove...\n";
            li_set.remove(ld_matrix[col]);
            // std::cout<<"\tRemoving: "<<ld_matrix[col];
            row_offset--;
        }
    }
    return li_set;
}
int *Set::dependencies_pos() const{
    int nbr_dep = dependencies();
    int* deps = new int[nbr_dep];
    //init array full of zeros
    for(int dep{0}; dep < nbr_dep; dep++){
        deps[dep] = 0;
    }
    int dep_index{0};

    // std::cout<<"In remove_dependencies()\n";
    Set li_set{*this};
    Matrix ld_matrix{li_set.convert_to_matrix()};
    Matrix li_matrix{li_set.convert_to_matrix().rref()};
    // std::cout<<"\tli_matrix: ";
    // std::cout<<li_matrix;

    int row_offset{0};

    for(int col{0}; col < li_matrix.cols(); col++){
        // std::cout<<"\tLooking through matrix... row: "<<col+row_offset<<" col: "<<col<<" col - offset: "<<col - row_offset<<" le: "<<li_matrix.matrix(col + row_offset, col)<<"\n"; 

        if(col + row_offset == li_matrix.rows() - 1){
            // std::cout<<"Removing remaining entries\n";
            for(int remaining_col{col + 1}; remaining_col < li_matrix.cols(); remaining_col++){
                deps[dep_index] = remaining_col;
                dep_index++;
            }
            break;
        }

        if(li_matrix.matrix(col + row_offset, col) != 1){
            // std::cout<<"\tAbout to remove...\n";
            deps[dep_index] = col;
            dep_index++;
            // std::cout<<"\tRemoving: "<<ld_matrix[col];
            row_offset--;
        }
    }

    return deps;

}

Set Set::col_space(Matrix const &m){
    Set col_sp{m};
    // std::cout<<col_sp<<m.rref();
    return col_sp.remove_dependencies();
}
Set Set::null_space(Matrix const &m){
    Set tmp{m};
    Set null_sp{};
    // std::cout<<"In null_space()\n\n";
    // std::cout<<"\tAbout to do rref\n";
    Matrix m_rref = m.rref();
    // std::cout<<"\n\tDone rref: \n";
    // std::cout<<m_rref;

    // std::cout<<"\tReceiving dependencies...  nbr: ";
    int dependencies = tmp.dependencies(); 
    int* dependency_locations = tmp.dependencies_pos();
    // std::cout<<dependencies<<"\n";

    // for(int dep{0}; dep < dependencies; dep++){
        // std::cout<<"\t\t"<<dependency_locations[dep]<<"\n";
    // }

    // std::cout<<"\tDependencies received.\n";
    // std::cout<<"\n\nDealing with dependencies...\n\n";

    for(int dep{0}; dep < dependencies; dep++){
        
        // std::cout<<"\t\tDealing with dep: "<<dep<<"...\n";
        Matrix m{m_rref.cols(), 1};
        // std::cout<<"\t\tCreated matrix m.\n";

        int dep_loc = dependency_locations[dep]; //this is the location of the first dependency
        // std::cout<<"\t\tGot dependency location: "<<dep_loc<<"\n";
        //start at the origin
        // std::cout<<"\t\tAbout to loop through matrix...\n";

        for(int row{0}, col{0}; col < m_rref.cols(); col++){

            // std::cout<<"\t\t\tOn col: "<<col<<" row: "<<row<<"\n";
            while(col != m_rref.get_le_col(row) && row != m_rref.rows() - 1){
                
                // std::cout<<"\t\t\t\tAssigning val... in while\n";

                if(col == dep_loc){ 
                    // std::cout<<"\t\t\t\tIn dep col: "<<col<<"\n";
                    m.matrix_[col][0] = 1;
                }
                else{
                    m.matrix_[col][0] = -m_rref.matrix(row, dep_loc);
                    if(m.matrix(col, 0) == -0) m.matrix_[col][0] = 0.0f; //handle the case where we have -0
                }

                col++;
                
            }
            // std::cout<<"\t\t\t\tAssigning val... after\n";
            if(col == dep_loc){ 
                // std::cout<<"\t\t\t\tIn dep col: "<<col<<"\n";
                m.matrix_[col][0] = 1;
            }
            else{
                m.matrix_[col][0] = -m_rref.matrix(row, dep_loc);
                if(m.matrix(col, 0) == -0) m.matrix_[col][0] = 0.0f; //handle the case where we have -0
            }

            // std::cout<<"\t\t\t\tAbout to inc row if needed\n";

            if(row < m_rref.rows() - 1) {
                // std::cout<<"\t\t\t\tInc row\n";
                row++;
            } //only increment rows if we haven't gotten to the last row
            // std::cout<<"Moving to next... m: "<<m;
        }

        // std::cout<<"\t\tDone parsing matrix. Inserting m into null_sp.\n";
        null_sp.insert(m);

    }

    // std::cout<<null_sp;


    //delete deps
    for(int dep{0}; dep < dependencies; dep++){
        dependency_locations[dep] = 0;
    }
    delete[] dependency_locations;
    dependency_locations = nullptr;

    return null_sp;
}

int Set::dim() const{
    Set tmp{remove_dependencies()};
    int dim{0};
    for(Node* p_ptr{p_head_}; p_ptr != nullptr; p_ptr = p_ptr->next(), dim++);
    return dim;
}   















