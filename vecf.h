//Vectors class with fixed size
#ifndef NL_VECF_H
#define NL_VECF_H

#include <iostream>
#include <limits>
#include <cmath>

using namespace std;

//Vector class
template <typename T>
class VecF{

public:

    //Data members
    int n; //number of values
    T* v; //values

    //Constructors, copy constructor, destructor
    VecF();
    VecF(int size);
    VecF(const VecF& source);
    ~VecF();

    //Member functions
    inline bool equals(const T& a, const T& b); //check for equality

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
    VecF operator + (const T& k);
    VecF operator - (const T& k);
    VecF operator * (const T& k);
    VecF operator / (const T& k);
    //Binary operators with vector
    void operator += (const VecF& source);
    void operator -= (const VecF& source);
    void operator *= (const VecF& source);
    void operator /= (const VecF& source);
    bool operator == (const VecF& source);
    VecF& operator = (const VecF& source);
    VecF operator + (const VecF& source);
    VecF operator - (const VecF& source);
    VecF operator * (const VecF& source);
    VecF operator / (const VecF& source);
    //Unary operators
    VecF operator - ();

    //Stream
    friend ostream& operator << (ostream &output, const VecF& source) {
        for (int i = 0; i < source.n; ++i) output << source.v[i] << endl;
        return output;
    };
};

#include "vecf.tpp"

#endif //NL_VECF_H
