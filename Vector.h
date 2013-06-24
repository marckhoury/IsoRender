#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

struct Point
{
	float x;
	float y;
	float z;
};

class Vector
{
private:
	float coord[3];
public:
	Vector(float x = 0, float y = 0, float z = 0);	
	Vector(const Vector& rhs);
	Vector operator+(Vector rhs);
	Vector operator-(Vector rhs);
	Vector operator-();
	Vector operator*(float c);
	float Dot(Vector& rhs);
	Vector Cross(Vector& rhs);
	void Normalize();
	float Length();
	float& operator[](int index);	
	void Swap(Vector& other);
	Vector& operator=(Vector rhs);
};
Vector operator*(float c, Vector rhs);

#endif
