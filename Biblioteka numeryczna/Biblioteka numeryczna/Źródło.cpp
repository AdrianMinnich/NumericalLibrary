#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include "Punkt.h"
#include <cstdlib>
//#include <cmath>
#include <iomanip>

using namespace std;

void wczytaj(vector<Punkt>& punkty, int pick);
double interpolacjaLagrangea(vector<Punkt> punkty, int x);
double interpolacjaNewtona(vector<Punkt> punkty, int x);
double ilorazRoznicowy(vector<Punkt> punkty, int lewy, int prawy);
void obliczWyznacznik();
void interpolacjaWielomianowa();
void calkowanieProstokatne();
double f(double x);
void calkowanieTrapezowe();
//void wczytaj(vector<double>& wspolczynniki, double &xp, double &xk, double &n);
//double f1(vector<double>w);
void calkowanieSimpson();
void calkowanieMonteCarlo();
void kwadraturaGaussa2D();
void metodaBisekcji();
void metodaNR();
double fp(double x);
void metodaRK1();
double f(double x, double y);
void metodaRK2();
void metodaRK4();


int main()
{
	vector<Punkt> punkty;
	int wybor;
	cout << "1. Interpolacja Lagrange'a\n2. Interpolacja Newtona\n3. Interpolacja wielomianowa (wzory Cramera)\n4. Calkowanie metoda prostokatow\n5. Calkowanie metoda trapezow\n6. Calkowanie metoda Simpsona\n7. Calkowanie metoda Monte Carlo\n8. Kwadratura Gaussa 2D\n9. Metoda bisekcji\n10. Metoda NR\n11. Metoda RK1\n12. Metoda RK2\n13. Metoda RK4\n0. Wyjscie\n\n";
	do {
		cout << "Wybor: ";
		cin >> wybor;
		switch (wybor)
		{
		case 1:
		{
			wczytaj(punkty, wybor);
			cout << "Podaj x: ";
			int x;
			cin >> x;
			cout << "Interpolacja Lagrange'a = " << interpolacjaLagrangea(punkty, x) << endl << endl;
			punkty.clear();
			break;
		}
		case 2:
		{
			wczytaj(punkty, wybor);
			cout << "Podaj x: ";
			int x;
			cin >> x;
			cout << "Interpolacja Newtona = " << interpolacjaNewtona(punkty, x) << endl << endl;
			punkty.clear();
			break;
		}
		case 3:
			interpolacjaWielomianowa();
			break;
		case 4:
			calkowanieProstokatne();
			break;
		case 5:
			calkowanieTrapezowe();
			break;
		case 6:
			calkowanieSimpson();
			break;
		case 7:
			calkowanieMonteCarlo();
			break;
		case 8:
			kwadraturaGaussa2D();
			break;
		case 9:
			metodaBisekcji();
			break;
		case 10:
			metodaNR();
			break;
		case 11:
			metodaRK1();
			break;
		case 12:
			metodaRK2();
			break;
		case 13:
			metodaRK4();
		default:
			break;
		}
	} while (wybor != 0);
	system("pause");
}

void wczytaj(vector<Punkt>& punkty, int pick)
{
	fstream plik;
	int x, y;
	if (pick == 1)
	{
		plik.open("Dane1.txt", ios::in);
		if (plik.good())
		{
			while (plik >> x >> y)
			{
				punkty.push_back(Punkt(x, y));
			}
		}
	}
	else if (pick == 2)
	{
		plik.open("Dane2.txt", ios::in);
		if (plik.good())
		{
			while (plik >> x >> y)
			{
				punkty.push_back(Punkt(x, y));
			}
		}
	}
}

double interpolacjaLagrangea (vector<Punkt> punkty, int x)
{
	double l;
	double L = 0.0;

	for (int i = 0; i < punkty.size(); i++)
	{
		l = 1.0;
		for (int j = 0; j < punkty.size(); j++)
		{
			if (j != i)
			{
				l *= ((x - punkty[j].x) / (punkty[i].x - punkty[j].x));
			}
		}
		L += l * punkty[i].y;
	}
	return L;
}

double interpolacjaNewtona(vector<Punkt> punkty, int x)
{
	double W = punkty[0].y;
	for (int i = 1; i < punkty.size(); i++)
	{
		double a = ilorazRoznicowy(punkty, 0, i); // iloraz roznicowy
		double m = 1;                             // (x-xn)
		for (int j = 0; j < i; j++)
		{
			m *= (x - punkty[j].x);
		}
		W += a * m;
	}
	return W;
}

double ilorazRoznicowy(vector<Punkt> punkty, int lewy, int prawy)
{
	double x1 = punkty[lewy].x;
	double x2 = punkty[prawy].x;
	double y1, y2;
	if (prawy - lewy == 1)
	{
		y1 = punkty[lewy].y;
		y2 = punkty[prawy].y;	
	}
	else
	{
		y1 = ilorazRoznicowy(punkty, lewy, prawy - 1);
		y2 = ilorazRoznicowy(punkty, lewy + 1, prawy);
	}
	return (y2 - y1) / (x2 - x1);
}

void obliczWyznacznik()
{
	double macierz[3][3];
	macierz[0][0] = 2;
	macierz[0][1] = 5;
	macierz[0][2] = 3;
	macierz[1][0] = 4;
	macierz[1][1] = 2;
	macierz[1][2] = 5;
	macierz[2][0] = 3;
	macierz[2][1] = 8;
	macierz[2][2] = 4;
	double wynik[3];
	wynik[0] = 5;
	wynik[1] = 4;
	wynik[2] = 9;

	double W = 0;
	double Wx = 0, Wy = 0, Wz = 0;
//	double suma = 0;
//	double iloczyn = 1;

	W = (macierz[0][0] * macierz[1][1] * macierz[2][2]) + (macierz[0][1] * macierz[1][2] * macierz[2][0]) + (macierz[0][2] * macierz[1][0] * macierz[2][1]) -
		(macierz[0][2] * macierz[1][1] * macierz[2][0]) - (macierz[0][0] * macierz[1][2] * macierz[2][1]) - (macierz[0][1] * macierz[1][0] * macierz[2][2]);

	Wx = wynik[0] * macierz[1][1] * macierz[2][2] + macierz[0][1] * macierz[1][2] * wynik[2] + macierz[0][2] * wynik[1] * macierz[2][1] -
		macierz[0][2] * macierz[1][1] * wynik[2] - wynik[0] * macierz[1][2] * macierz[2][1] - macierz[0][1] * wynik[1] * macierz[2][2];

	Wy = macierz[0][0] * wynik[1] * macierz[2][2] + wynik[0] * macierz[1][2] * macierz[2][0] + macierz[0][2] * macierz[1][0] * wynik[2] -
		macierz[0][2] * wynik[1] * macierz[2][0] - macierz[0][0] * macierz[1][2] * wynik[2] - wynik[0] * macierz[1][0] * macierz[2][2];

	Wz = macierz[0][0] * macierz[1][1] * wynik[2] + macierz[0][1] * wynik[1] * macierz[2][0] + wynik[0] * macierz[1][0] * macierz[2][1] -
		wynik[0] * macierz[1][1] * macierz[2][0] - macierz[0][0] *wynik[1] * macierz[2][1] - macierz[0][1] * macierz[1][0] * wynik[2];


	double x = Wx / W;
	double y = Wy / W;
	double z = Wz / W;

	cout << "x = " << x << endl << "y= " << y << endl << "z= " << z << endl;
}

void interpolacjaWielomianowa()
{
	double macierz[3][3];
	/*
	macierz[0][0] = pow(2, 2);
	macierz[0][1] = pow(2, 1);
	macierz[0][2] = pow(2, 0);
	macierz[1][0] = pow(3, 2);
	macierz[1][1] = pow(3, 1);
	macierz[1][2] = pow(3, 0);
	macierz[2][0] = pow(5, 2);
	macierz[2][1] = pow(5, 1);
	macierz[2][2] = pow(5, 0);
	*/
	macierz[0][0] = 4;
	macierz[0][1] = 2;
	macierz[0][2] = 1;
	macierz[1][0] = 16;
	macierz[1][1] = 4;
	macierz[1][2] = 1;
	macierz[2][0] = 36;
	macierz[2][1] = 6;
	macierz[2][2] = 1;
	double wynik[3];
	wynik[0] = 8;
	wynik[1] = 5;
	wynik[2] = 3;

	double W = 0;
	double Wx = 0, Wy = 0, Wz = 0;
	

	W = (macierz[0][0] * macierz[1][1] * macierz[2][2]) + (macierz[0][1] * macierz[1][2] * macierz[2][0]) + (macierz[0][2] * macierz[1][0] * macierz[2][1]) -
		(macierz[0][2] * macierz[1][1] * macierz[2][0]) - (macierz[0][0] * macierz[1][2] * macierz[2][1]) - (macierz[0][1] * macierz[1][0] * macierz[2][2]);

	Wx = wynik[0] * macierz[1][1] * macierz[2][2] + macierz[0][1] * macierz[1][2] * wynik[2] + macierz[0][2] * wynik[1] * macierz[2][1] -
		macierz[0][2] * macierz[1][1] * wynik[2] - wynik[0] * macierz[1][2] * macierz[2][1] - macierz[0][1] * wynik[1] * macierz[2][2];

	Wy = macierz[0][0] * wynik[1] * macierz[2][2] + wynik[0] * macierz[1][2] * macierz[2][0] + macierz[0][2] * macierz[1][0] * wynik[2] -
		macierz[0][2] * wynik[1] * macierz[2][0] - macierz[0][0] * macierz[1][2] * wynik[2] - wynik[0] * macierz[1][0] * macierz[2][2];

	Wz = macierz[0][0] * macierz[1][1] * wynik[2] + macierz[0][1] * wynik[1] * macierz[2][0] + wynik[0] * macierz[1][0] * macierz[2][1] -
		wynik[0] * macierz[1][1] * macierz[2][0] - macierz[0][0] * wynik[1] * macierz[2][1] - macierz[0][1] * macierz[1][0] * wynik[2];


	double a2 = Wx / W;
	double a1 = Wy / W;
	double a0 = Wz / W;

	//cout << "a2 = " << a2 << endl << "a1 = " << a1 << endl << "a0 = " << a0 << endl;
	cout << "Wielomian: " << a2 << "x^2 + " << a1 << "x + " << a0 << endl << endl;
	
}

void calkowanieProstokatne()
{
	double xp = 1, xk = 5, n = 4;
	double dx, s = 0;
	dx = (xk - xp) / n;
	for (int i = 1; i <= n; i++)
		s += f(xp + i);
	s *= dx;
	cout << "Wynik calkowania: " << s << endl << endl;
}

double f(double x)
{
	return x * x * x - 2;
}

double fp(double x)
{
	return 3 * x*x;
}

void calkowanieTrapezowe()
{
	double xp = 1, xk = 5, n = 4;
	double dx, s = 0;
	dx = (xk - xp) / n;
	for (int i = 1; i < n; i++)
		s += f(xp + i * dx);
	s = (s + (f(xp) + f(xk)) / 2) * dx;
	cout << "Wynik calkowania: " << s << endl << endl;
}

/*void wczytaj(vector<double>& wspolczynniki, double &xp, double &xk, double &n)
{
	string temp;
	fstream plik;
	double x;
	plik.open("Dane3.txt", ios::in);
	if (plik.good())
	{
		while (plik >> x)
			wspolczynniki.push_back(x);
	}
	n = wspolczynniki[wspolczynniki.size() - 1];
	xk = wspolczynniki[wspolczynniki.size() - 2];
	xp = wspolczynniki[wspolczynniki.size() - 3];

	wspolczynniki.pop_back;
	wspolczynniki.pop_back;
	wspolczynniki.pop_back;

}*/

/*double f1(vector<double>w)
{
	int stopien = w.size() - 2;
	while ()
	{

	}
}*/

void calkowanieSimpson()
{
	double s = 0, st = 0, dx, a = 0, x;
	double xp = 1, xk = 5, n = 4, h = 1;
	dx = (xk - xp) / n;
	for (int i = 1; i <= n; i = i++)
	{	
		x = xp + i * dx;
		st += f(x - dx / 2);
		if (i < n) s += f(x);
	}
	s = dx / 6 * (f(xp) + f(xk) + 2 * s + 4 * st);
	cout << "Wynik calkowania: " << s << endl << endl;
}

void calkowanieMonteCarlo()
{
	double s = 0, dx = 0;
	double x[4] = { 1.5, 2.6, 3.8, 4.5 };
	double xp = 1, xk = 4, n = 3;
	dx = abs(xk - xp);
	for (int i = 1; i <= n; i++)
		s += f(x[i - 1]) / n;
	s = s * dx;
	cout << "Wynik calkowania: " << s << endl << endl;
}

void kwadraturaGaussa2D()
{
	double wsp_x[4], wsp_y[4];
	double waga[2], punkt[2];
	double fksztalt[2][2][4];
	double poch_ksi[2][4], poch_ni[2][4], fun_detj[2][2];

	fstream plik;
	plik.open("input.txt", ios::in);
	if (plik.good())
	{
		/*for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{

			}
		}*/
		plik >> wsp_x[0] >> wsp_y[0] >> wsp_x[1] >> wsp_y[1] >> wsp_x[2] >> wsp_y[2] >> wsp_x[3] >> wsp_y[3];
	}
	waga[0] = 1.0;
	waga[1] = 1.0;
	punkt[0] = -0.8543743;
	punkt[1] = 0.8543743;

	for (int j = 0; j <= 1; j++)
	{
		for (int i = 0; i <= 1; i++)
		{
			fksztalt[i][j][0] = 0.25 * (1.0 - punkt[i]) * (1.0 - punkt[j]);
			fksztalt[i][j][1] = 0.25 * (1.0 + punkt[i]) * (1.0 - punkt[j]);
			fksztalt[i][j][2] = 0.25 * (1.0 + punkt[i]) * (1.0 + punkt[j]);
			fksztalt[i][j][3] = 0.25 * (1.0 - punkt[i]) * (1.0 + punkt[j]);

			poch_ksi[j][0] = -0.25 * (1.0 - punkt[j]);
			poch_ksi[j][1] = 0.25 * (1.0 - punkt[j]);
			poch_ksi[j][2] = 0.25 * (1.0 + punkt[j]);
			poch_ksi[j][3] = -0.25 * (1.0 + punkt[j]);
			poch_ni[i][0] = -0.25 * (1.0 - punkt[i]);
			poch_ni[i][1] = -0.25 * (1.0 + punkt[i]);
			poch_ni[i][2] = 0.25 * (1.0 + punkt[i]);
			poch_ni[i][3] = 0.25 * (1.0 - punkt[i]);
		}
	}

	double dxdksi, dydksi, dxdni, dydni;
	for (int j = 0; j <= 1; j++)
	{
		for (int i = 0; i <= 1; i++)
		{
			dxdksi = poch_ksi[j][0] * wsp_x[0] + poch_ksi[j][1] * wsp_x[1] + poch_ksi[j][2] * wsp_x[2] + poch_ksi[j][3] * wsp_x[3];
			dydksi = poch_ksi[j][0] * wsp_y[0] + poch_ksi[j][1] * wsp_y[1] + poch_ksi[j][2] * wsp_y[2] + poch_ksi[j][3] * wsp_y[3];
			dxdni = poch_ni[i][0] * wsp_x[0] + poch_ni[i][1] * wsp_x[1] + poch_ni[i][2] * wsp_x[2] + poch_ni[i][3] * wsp_x[3];
			dydni = poch_ni[i][0] * wsp_y[0] + poch_ni[i][1] * wsp_y[1] + poch_ni[i][2] * wsp_y[2] + poch_ni[i][3] * wsp_y[3];

			fun_detj[i][j] = dxdksi * dydni - dxdni * dydksi;
		}
	}

	double powierzchnia = 0.0 ;
	for (int j = 0; j <= 1; j++)
	{
		for (int i = 0; i <= 1; i++)
		{
			powierzchnia += powierzchnia + abs(fun_detj[i][j]) * waga[i] * waga[j];
		}
	}

	cout << "Powierzchnia wynosi = " << powierzchnia << endl << endl;

}

void metodaBisekcji()
{
	double a, b, x1, fa, fb, fx;
	cout << "a = ";
	cin >> a;
	cout << "b = ";
	cin >> b;
	fa = f(a);
	fb = f(b);
	if (fa * fb > 0)
		cout << "Sa tych samych znakow.\n";
	else
	{
		while (fabs(a - b) > 0.005)
		{
			x1 = (a + b) / 2;
			fx = f(x1);
			if (fa *fx < 0)
				b = x1;
			else
			{
				a = x1;
				fa = fx;
			}
		}
		cout << "x0 = " <<  x1 << endl << endl;
	}
}

void metodaNR()
{
	double x0, x1, f0, f1;
	int i;
	cout << "x0 = ";
	cin >> x0;
	x1 = x0 - 1;
	f0 = f(x0);
	i = 64;
	while (i && (fabs(x1 - x0) > 0.005) && (fabs(x0)) > 0.005)
	{
		f1 = fp(x0);
		if (fabs(f1) < 0.005)
		{
			cout << "Zly punkt\n";
			i = 0;
			break;
		}
		x1 = x0;
		x0 = x0 - f0 / f1;
		f0 = f(x0);
		if (!(--i))
			cout << "Przekroczono limit obiegów\n";
	}
	if (i)
		cout << "x0 = " << x0 << endl << endl;
}

void metodaRK1()
{
	double x0, y0, b, h, n;
	x0 = 0;
	y0 = 1;
	b = 3;
	h = 0.1;
	n = (b - x0) / h;
	double x1, y1;
	fstream plik;
	ofstream file;
	file.open("RK1h01e.csv");
	/*if(h == 0.01) plik.open("RK1h001.txt", ios::out | ios::app);
	else if(h == 0.1) plik.open("RK1h01.txt", ios::out | ios::app);
	else if(h == 0.5) plik.open("RK1h05.txt", ios::out | ios::app);*/
		
	for (int i = 0; i < n; i++)
	{
		y1 = y0 + h * f(x0, y0);
		x1 = x0 + h;	
		if (file.good())
		{
			//plik << setw(10) << left << x1 << setw(10) << left << y1 << endl;
			//plik << x1 << ";" << y1 << endl;
			file  << x0 << "," << y0 << "\n";
		}
		y0 = y1;
		x0 = x1;
		y1 = 0;
		x1 = 0; 
	}
	//plik.close();
	file.close();
	cout << "Pomyslnie zapisano do pliku.\n";
}

double f(double x, double y)
{
	return x+y;
}

void metodaRK2()
{
	double x0, y0, b, h, n;
	x0 = 0;
	y0 = 2;
	b = 0.3;
	h = 0.02;
	n = (b - x0) / h;
	double x1, y1;
	ofstream file;
	file.open("RK2h02.csv");
	for (int i = 0; i < n; i++)
	{
		y1 = y0 + h / 2 * (f(x0, y0) + f(x0 + h, y0 + h * f(x0, y0)));
		x1 = x0 + h;
		if (file.good())
		{
			file << x0 << "," << y0 << "\n";
		}
		y0 = y1;
		x0 = x1;
		y1 = 0;
		x1 = 0;
	}
	file.close();
	cout << "Pomyslnie zapisano do pliku.\n";
}

void metodaRK4()
{
	
	double x0, y0, b, h, n, fi;
	x0 = 0;
	y0 = 1;
	b = 0.2;
	h = 0.06;
	n = (b - x0) / h;
	double x1, y1;
	ofstream file;
	double k1, k2, k3, k4;
	file.open("RK4h006.csv");
	for (int i = 1; i < n; i++)
	{
		k1 = f(x0, y0);
		k2 = f((x0 + 0.5*h), (y0 + 0.5*h*k1));
		k3 = f((x0 + 0.5*h), (y0 + 0.5*h*k2));
		k4 = f((x0 + h), (y0 + h * k3));
		fi = (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4);
		if (file.good())
		{
			file << x0 << "," << y0 << ", " << k1 << "\n";
			file << x0 + 0.5 * h << "," << y0 + 0.5 * k1 << ", " << k2 << "\n";
			file << x0 + 0.5 * h << "," << y0 + 0.5 * k2 << ", " << k3 << "\n";
			file << x0 + h << "," << y0 + k3 << ", " << k4 << "\n";
		}

		y1 = y0 + (1 / 6) * (k1 + 2*k2 + 2*k3 + k4);
		x1 = x0 + h;

		y0 = y1;
		x0 = x1;
		y1 = 0;
		x1 = 0;

	}
	file.close();
	cout << "Pomyslnie zapisano do pliku.\n";


}
