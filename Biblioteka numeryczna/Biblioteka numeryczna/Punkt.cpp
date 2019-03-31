#include "Punkt.h"



Punkt::Punkt()
{
}

Punkt::Punkt(double xx, double yy)
{
	x = xx;
	y = yy;
}


Punkt::~Punkt()
{
}

void Punkt::set(double xx, double yy)
{
	x = xx;
	y = yy;
}

double Punkt::getX()
{
	return x;
}

double Punkt::getY()
{
	return y;
}
