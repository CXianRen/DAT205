#pragma once

class Vec3
{
public:
  float n[3];

  // Constructors
  Vec3();
  Vec3(const float x, const float y, const float z);
  Vec3(const Vec3 &v); // copy constructor

  // Assignment operators
  Vec3 &operator=(const Vec3 &v);   // assignment of a Vec3
  Vec3 &operator+=(const Vec3 &v);  // incrementation by a Vec3
  Vec3 &operator-=(const Vec3 &v);  // decrementation by a Vec3
  Vec3 &operator*=(const float d); // multiplication by a constant
  Vec3 &operator/=(const float d); // division by a constant
  float &operator[](int i);        // indexing
  float operator[](int i) const;   // read-only indexing
  void set(const float x, const float y, const float z);

  // special functions
  float norm() const;             // length of a Vec3
  float sqrLength() const;        // squared length of a Vec3
  Vec3 &normalize();               // normalize a Vec3 in place
  Vec3 cross(const Vec3 &v) const; // cross product: self cross v

  // friends
  friend Vec3 operator-(const Vec3 &v);                    // -v1
  friend Vec3 operator+(const Vec3 &a, const Vec3 &b);     // v1 + v2
  friend Vec3 operator-(const Vec3 &a, const Vec3 &b);     // v1 - v2
  friend Vec3 operator*(const Vec3 &a, const float d);    // v1 * scalar
  friend Vec3 operator*(const float d, const Vec3 &a);    // scalar * v1
  friend Vec3 operator*(const Vec3 &a, const Vec3 &b);     // piecewise muliply
  friend Vec3 operator/(const Vec3 &a, const float d);    // v1 / scalar
  friend Vec3 operator^(const Vec3 &a, const Vec3 &b);     // cross product
  friend bool operator==(const Vec3 &a, const Vec3 &b);    // v1 == v2 ?
  friend bool operator!=(const Vec3 &a, const Vec3 &b);    // v1 != v2 ?
  friend Vec3 prod(const Vec3 &a, const Vec3 &b);          // term by term *
  friend float dot(const Vec3 &a, const Vec3 &b);         // dot product
  friend float distance(const Vec3 &a, const Vec3 &b);    // distance
  friend float distanceSqr(const Vec3 &a, const Vec3 &b); // distance sqr
};