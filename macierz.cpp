#include <iostream>
#include <fstream>

using namespace std;
ofstream plik;

double** TworzMacierz(int IloscWierszy, int IloscKolumn)
{

    double **Tablica = new double *[IloscWierszy];
    for(int i = 0; i < IloscWierszy ; i++)
    {
        Tablica[i] = new double [IloscKolumn];
    }

    // Zerowanie nowo utworzonej tablicy
    for(int i = 0; i < IloscWierszy ; i++)
    {
        for(int k = 0; k < IloscKolumn ; k++)
        {
            Tablica[i][k] = 0;
        }
    }

    return Tablica;
}

void DrukujMacierz(double **Tablica, int IloscWierszy, int IloscKolumn, double x, double h, double d, double dt)
{

    plik.open("MojaTablica.csv");
    plik << "\t";
    for(int k = 0; k < IloscKolumn ; k++)
    {
        plik << x << ";";
        x +=h;
    }plik << endl;
    for(int i = 0; i < IloscWierszy ; i++)
    {
        plik << d << ";";
        for(int k = 0; k < IloscKolumn ; k++)
        {
            // cout << Tablica[i][k] << " ";
            plik << Tablica[i][k] << ";";
        }plik << endl;
        d += dt;
    }

    plik.close();
}