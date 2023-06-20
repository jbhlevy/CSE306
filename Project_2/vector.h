//
//  vector.h
//  Project_1
//
//  Created by John Levy on 07/05/2023.
//

#ifndef vector_h
#define vector_h

#include <vector>
#include <cmath>

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
    void operator+=(const Vector& other)
    {
        data[0] += other[0];
        data[1] += other[1];
        data[2] += other[2];
    }
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector operator*(const Vector& a, const Vector& b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] ; // ; + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector operator-(const Vector& a){
    return Vector(-a[0], -a[1], -a[2]);
}

Vector findBarycenter(const Vector& A, const Vector& B, const Vector& C){
    Vector P = (A + B + C )/3;
    return P;
}

//Not a vector but useful here to avoid unecessary redundant inclusions of file
double sqr(double x){
    return std::pow(x, 2);
}

#endif /* vector_h */
