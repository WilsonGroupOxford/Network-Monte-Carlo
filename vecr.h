//Vector class with reserved size
#ifndef NL_VECR_H
#define NL_VECR_H

#include <iostream>
#include <limits>
#include <cmath>

using namespace std;

/* Vector class with reserved size
 * Maximum size defines limit on the number of values
 * Current size determines accessible number of values */
template <typename T>
class VecR{

public:

    //Data members
    int n, nMax; //number of values, maximum number of values
    T* v; //values

    //Constructors, copy constructor, destructor
    VecR();
    VecR(int maxSize);
    VecR(int size, int maxSize);
    VecR(const VecR& source);
    ~VecR();

    //Member functions
    inline bool equals(const T& a, const T& b); //check for equality
    void setSize(int size); //set current size
    void resetMaxSize(int maxSize); //reset maximum size
    void addValue(T value); //add value to end
    void delValue(T value); //remove first instance of value from vector
    void swapValue(T vDel, T vAdd, bool swapAll=true); //swap a value in place of another
    void swapValue(T vDel, T vAdd, T vBetween0, T vBetween1); //swap value in place of another at specific point
    void insertValue(T vInsert, T vBetween0, T vBetween1); //insert a value between two others

    //Subscript operator
    T& operator [] (int i);
    //Binary operators with constant
    void operator = (const T& k);
    void operator += (const T& k);
    void operator -= (const T& k);
    void operator *= (const T& k);
    void operator /= (const T& k);
    bool operator == (const T& k);
    bool operator < (const T& k);
    bool operator > (const T& k);
    VecR operator + (const T& k);
    VecR operator - (const T& k);
    VecR operator * (const T& k);
    VecR operator / (const T& k);
    //Binary operators with vector
    void operator += (const VecR& source);
    void operator -= (const VecR& source);
    void operator *= (const VecR& source);
    void operator /= (const VecR& source);
    bool operator == (const VecR& source);
    VecR& operator = (const VecR& source);
    VecR operator + (const VecR& source);
    VecR operator - (const VecR& source);
    VecR operator * (const VecR& source);
    VecR operator / (const VecR& source);
    //Unary operators
    VecR operator - ();

    //Stream
    friend ostream& operator << (ostream &output, const VecR& source) {
        for (int i = 0; i < source.n; ++i) output << source.v[i] << endl;
        return output;
    };
};

#include "vecr.tpp"

#endif //NL_VECR_H
