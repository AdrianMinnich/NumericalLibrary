#pragma once
class Punkt
{
public:
	double x;
	double y;
public:
	Punkt();
	Punkt(double xx, double yy);
	~Punkt();
	void set(double xx, double yy);
	double getX();
	double getY();
	
};

