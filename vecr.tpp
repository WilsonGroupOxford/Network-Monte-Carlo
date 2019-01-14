#include "vecr.h"

//##### VECTOR RESERVED SIZE CLASS #####

//Default constructor
template <typename T>
VecR<T>::VecR() {
    this->n = 0;
    this->nMax = 1;
    this->v = new T[this->nMax]();
}

//Construct with maximum size
template <typename T>
VecR<T>::VecR(int maxSize) {
    if(maxSize<0) throw "Vector cannot be instantiated with negative size";
    this->n = maxSize;
    this->nMax = maxSize;
    this->v = new T[this->nMax]();
}

//Construct with size and maximum size
template <typename T>
VecR<T>::VecR(int size, int maxSize) {
    if(maxSize<0) throw "Vector cannot be instantiated with negative size";
    this->n = size;
    this->nMax = maxSize;
    this->v = new T[this->nMax]();
}

//Construct from VecR
template <typename T>
VecR<T>::VecR(const VecR& source){
    this->n = source.n;
    this->nMax = source.nMax;
    this->v = new T[this->nMax]();
    for (int i=0; i<this->n; ++i) this->v[i] = source.v[i];
}

//Destructor, clear allocated memory
template <typename T>
VecR<T>::~VecR() {
    delete[] this->v;
}

//Check for equality of any value type
template <typename T>
bool VecR<T>::equals(const T& a, const T& b) {
    return a==b || abs(a-b)<abs(min(a,b))*numeric_limits<T>::epsilon();
}

//Set current size of vector
template <typename T>
void VecR<T>::setSize(int size) {
    if(size>this->nMax) throw "Cannot set vector size larger than reserved size";
    this->n = size;
}

//Reset maximum size of vector
template <typename T>
void VecR<T>::resetMaxSize(int maxSize) {
    if(maxSize!=this->nMax){
        T* vCopy = new T[this->n]();
        for(int i=0; i<this->n; ++i) vCopy[i]=this->v[i];
        delete[] this->v;
        this->nMax = maxSize;
        this->v=new T[this->nMax]();
        for(int i=0; i<this->n; ++i) this->v[i]=vCopy[i];
        delete[] vCopy;
    }
}

//Add value to end of vector
template <typename T>
void VecR<T>::addValue(T value) {
    if(this->n == this->nMax) throw "Cannot add to vector to make larger than reserved size";
    this->v[this->n] = value;
    ++this->n;
}

//Deletes first instance of value from vector
template <typename T>
void VecR<T>::delValue(T value) {
    bool del=false;
    int d;
    for(int i=0; i<this->n; ++i){
        if(equals(this->v[i],value)){
            d=i;
            del=true;
            break;
        }
    }
    if(!del) throw "Cannot delete value as not present in vector";
    for(int i=d; i<this->n-1; ++i) this->v[i]=this->v[i+1];
    --this->n;
}

//Swaps value from vector
template <typename T>
void VecR<T>::swapValue(T vDel, T vAdd, bool swapAll) {
    bool swap=false;
    if(swapAll){
        for(int i=0; i<this->n; ++i){
            if(equals(this->v[i],vDel)){
                this->v[i]=vAdd;
                swap=true;
            }
        }
    }
    else{
        for(int i=0; i<this->n; ++i){
            if(equals(this->v[i],vDel)){
                this->v[i]=vAdd;
                swap=true;
                break;
            }
        }
    }
    if(!swap) throw "Cannot swap value as not present in vector";
}

//Insert value in vector between two others
template <typename T>
void VecR<T>::insertValue(T vInsert, T vBetween0, T vBetween1) {
    if(this->n == this->nMax) throw "Cannot insert to vector to make larger than reserved size";
    bool insert=false;
    int insertPos=-1;
    for(int i=0; i<this->n; ++i){
        if(this->v[i]==vBetween0 && this->v[(i+1)%this->n]==vBetween1){
            insertPos=(i+1)%this->n;
            insert=true;
        }
        else if(this->v[i]==vBetween1 && this->v[(i+1)%this->n]==vBetween0){
            insertPos=(i+1)%this->n;
            insert=true;
        }
        if(insert) break;
    }
    if(!insert) throw "Cannot insert value as surrounding values not present in vector";
    for(int i=this->n; i>insertPos; --i) this->v[i]=this->v[i-1];
    this->v[insertPos]=vInsert;
    ++this->n;
}

//Subscript operator
template <typename T>
T& VecR<T>::operator[](int i) {
    if(this->n<=i) throw "Vector subscript out of bounds";
    return this->v[i];
}

//Binary Operators with constant
template <typename T>
void VecR<T>::operator=(const T& k) {
    for(int i=0; i<this->n; ++i) this->v[i] = k;
}

template <typename T>
void VecR<T>::operator+=(const T& k) {
    for(int i=0; i<this->n; ++i) this->v[i] += k;
}

template <typename T>
void VecR<T>::operator-=(const T& k) {
    for(int i=0; i<this->n; ++i) this->v[i] -= k;
}

template <typename T>
void VecR<T>::operator*=(const T& k) {
    for(int i=0; i<this->n; ++i) this->v[i] *= k;
}

template <typename T>
void VecR<T>::operator/=(const T& k) {
    for(int i=0; i<this->n; ++i) this->v[i] /= k;
}

template <typename T>
bool VecR<T>::operator==(const T& k) {
    bool equality=true;
    for(int i=0; i<this->n; ++i){
        if(!equals(this->v[i],k)) equality=false;
    }
    return equality;
}

template <typename T>
bool VecR<T>::operator<(const T& k) {
    bool lt=true;
    for(int i=0; i<this->n; ++i){
        if(this->v[i]>=k) lt=false;
    }
    return lt;
}

template <typename T>
bool VecR<T>::operator>(const T& k) {
    bool gt=true;
    for(int i=0; i<this->n; ++i){
        if(this->v[i]<=k) gt=false;
    }
    return gt;
}

template <typename T>
VecR<T> VecR<T>::operator+(const T& k) {
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]+k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator-(const T& k) {
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]-k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator*(const T& k) {
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]*k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator/(const T& k) {
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]/k;
    return vec;
}

//Binary operators with VecR
template <typename T>
void VecR<T>::operator+=(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    for (int i = 0; i < this->n; ++i) this->v[i] += source.v[i];
}

template <typename T>
void VecR<T>::operator-=(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    for (int i = 0; i < this->n; ++i) this->v[i] -= source.v[i];
}

template <typename T>
void VecR<T>::operator*=(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    for (int i = 0; i < this->n; ++i) this->v[i] *= source.v[i];
}

template <typename T>
void VecR<T>::operator/=(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    for (int i = 0; i < this->n; ++i) this->v[i] /= source.v[i];
}

template <typename T>
bool VecR<T>::operator==(const VecR& source) {
    bool equality=true;
    if(this->n!=source.n) equality=false;
    else{
        for(int i=0; i<this->n; ++i){
            if(!equals(this->v[i],source.v[i])) equality=false;
        }
    }
    return equality;
}

template <typename T>
VecR<T>& VecR<T>::operator=(const VecR& source) {
    if(this != &source){
        delete[] this->v;
        this->n=source.n;
        this->nMax=source.nMax;
        this->v=new T[this->nMax]();
        for(int i=0; i<this->n; ++i) this->v[i]=source.v[i];
    }
    return *this;
}

template <typename T>
VecR<T> VecR<T>::operator+(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]+source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator-(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]-source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator*(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]*source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator/(const VecR& source) {
    if(this->n != source.n) throw "Cannot perform binary operation on vectors of different sizes";
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=this->v[i]/source.v[i];
    return vec;
}

//Unary operators
template <typename T>
VecR<T> VecR<T>::operator-(){
    VecR<T> vec(this->n);
    for(int i=0; i<this->n; ++i) vec[i]=-this->v[i];
    return vec;
}
