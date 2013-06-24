#include "Vector.h"

Vector::Vector(float x, float y, float z)
{
	coord[0] = x;
	coord[1] = y;
	coord[2] = z;
}	
	
Vector::Vector(const Vector& rhs)
{
	coord[0] = rhs.coord[0];
	coord[1] = rhs.coord[1];
	coord[2] = rhs.coord[2];
}

Vector Vector::operator+(Vector rhs)
{
	Vector lhs;
	lhs.coord[0] = coord[0] + rhs.coord[0];
	lhs.coord[1] = coord[1] + rhs.coord[1];
	lhs.coord[2] = coord[2] + rhs.coord[2];
	return lhs;
}
	
Vector Vector::operator-(Vector rhs)
{
	Vector lhs;
	lhs.coord[0] = coord[0] - rhs.coord[0];
	lhs.coord[1] = coord[1] - rhs.coord[1];
	lhs.coord[2] = coord[2] - rhs.coord[2];
	return lhs;
}	

Vector Vector::operator-()
{
	Vector temp;
	temp.coord[0] = -coord[0];
	temp.coord[1] = -coord[1];
	temp.coord[2] = -coord[2];
	return temp;
}

Vector Vector::operator*(float c)
{
	Vector temp;
	temp[0] = c * coord[0];
	temp[1] = c * coord[1];
	temp[2] = c * coord[2];
	return temp;
}	

float Vector::Dot(Vector& rhs)
{
	return coord[0] * rhs.coord[0] + coord[1] * rhs.coord[1] + coord[2] * rhs.coord[2];
}

Vector Vector::Cross(Vector& rhs)
{
	Vector temp;
	temp.coord[0] = (coord[1]*rhs.coord[2] - coord[2]*rhs.coord[1]);
	temp.coord[1] = (coord[2]*rhs.coord[0] - coord[0]*rhs.coord[2]);
	temp.coord[2] = (coord[0]*rhs.coord[1] - coord[1]*rhs.coord[0]);
	return temp;
}
	
void Vector::Normalize()
{
	float length = sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2]);
	if(length == 0)
	{return;}
	coord[0] = coord[0]/length;
	coord[1] = coord[1]/length;
	coord[2] = coord[2]/length;
}

float Vector::Length()
{
	return  sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2]);
}	

float& Vector::operator[](int index)
{
	return coord[index];
}	
	
void Vector::Swap(Vector& other)
{
	float temp;
	temp = coord[0];
	coord[0] = other.coord[0];
	other.coord[0] = temp;
	
	temp = coord[1];
	coord[1] = other.coord[1];
	other.coord[1] = temp;
	
	temp = coord[2];
	coord[2] = other.coord[2];
	other.coord[2] = temp;
}

Vector& Vector::operator=(Vector rhs)
{
	if(this != &rhs)
	{
		coord[0] = rhs.coord[0];
		coord[1] = rhs.coord[1];
		coord[2] = rhs.coord[2];
	}
	return *this;
}

Vector operator*(float c, Vector rhs)
{
	Vector temp;
	temp[0] = c * rhs[0];
	temp[1] = c * rhs[1];
	temp[2] = c * rhs[2];
	return temp;
}
