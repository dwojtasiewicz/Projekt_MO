#include "macierz.h"
#include "calerf.h"
#include <cstdio>
#include <fstream>
#include <sys/timeb.h>
#include <iostream>
#include <math.h>
using namespace std;

/*zmienne plikowe*/
ofstream Analitycznie;
ofstream Pomiary;
ofstream MaxBledy;

/*zmienne uzyte w programie:*/

/*parametr podany*/
double Lambda;
/*poczatkowe t*/
double t0;
/*maksymalne t*/
double tMax;
/*zmienna przydatna przy ustalaniu przedzialu*/
double a;
/*parametr podany*/
double D;
/*krok*/
double h;
/*parametr podany*/
double b;
/*przyrost czasowy*/
double dt;
/*liczba wierszy (double)*/
double d;
/*nasza tablica wynikowa*/
double **Tablica;
/*zmienne przetrzymujace ilosc kolumn i wierszy*/
int IloscKolumn;
int IloscWierszy;


/*metoda inicjujaca dane dla N kolumn*/
void inicjujDane(int N)
{
    /*ustalamy podana ilosc kolumn*/
    IloscKolumn= N;
    /*ustalamy parametr lambda (dla metod bezposrednich jest rowny 0.4
    - w naszym przypadku jest rowny 1*/
    Lambda= 0.4;
    /*ustalamy t poczatkowe*/
    t0= 0.0;
    /*ustalamy t koncowe*/
    tMax=2.0;
    /*ustalamy parametry b i D*/
   // b= 0.1;
    D= 1.0;

    /*obliczamy a, h, t, dt oraz liczbe wierszy macierzy*/
    a= 6* sqrt( D* tMax );
    h= (a+ a)/ (IloscKolumn- 1);
    dt= (Lambda* h* h)/ D;
    d= (tMax/ dt)+ 1;
    IloscWierszy= (int)d;
}

/*metoda liczaca wektor Ni dla algorytmu Thomasa*/
double *TworzNi(double *l, double *d, double *u)
{
    double *ni= new double[IloscKolumn];
    ni[0]= d[0];
    for(int i= 1; i< IloscKolumn; i++)
        ni[i]= d[i]- ( l[i]* u[i- 1] )/ ni[i- 1];
    return ni;
}

/*metoda Thomasa*/
void Thomas(double *l, double *ni, double *u, double*x, double *B)
{
    /*wektor Ni zostal juz obliczony*/


    /*teraz obliczamy wektor r*/
    double *r= new double[IloscKolumn];
    r[0]= B[0];
    for(int i= 1; i< IloscKolumn; i++)
        r[i]= B[i]- ( l[i]* r[i- 1] )/ ni[i- 1];

    /*na koniec zostalo nam wyznaczyc wektor wynikowy X*/
    x[IloscKolumn- 1]= r[IloscKolumn- 1]/ ni[IloscKolumn- 1];
    for(int i= IloscKolumn- 2; i>= 0; i--)
        x[i] = ( r[i]- u[i]* x[i+ 1] )/ ni[i];
    delete[] r;
}

/*rownanie analityczne*/
double RownanieAnalityczne(double x, double t)
{

    return erfc(x/(2*sqrt(D*t)));
}

/*metoda ustalajaca warunki poczatkowe i brzegowe*/
void Warunki()
{
    /*poczatek przedzialu x*/
    double x= 0;
    /*obliczamy przyrost*/
    double Delta_x= (a)/ (IloscKolumn- 1);
    /*ustalamy warunki poczatkowe*/
    for(int i= 0; i< IloscKolumn; i++)
    {
            Tablica[0][i]= 0;
    }
    /*ustalamy warunki brzegowe*/
    for(int j= 0; j< IloscWierszy; j++)
    {
        Tablica[j][0]= 1;
        Tablica[j][IloscKolumn-1]= 0;
    }
}

/*metoda wypisujaca parametry dla aktualnego dzialania programu*/
void WypiszDane()
{
    cout << "_____________________________________________________________________________"
         << "\n\nParametry zadania:\n"
         << "Lambda: "<< Lambda<< "\n"
         << "t poczatkowe: "<< t0<< "\n"
         << "t koncowe: "<< tMax<< "\n"
         << "a: "<< a<< "\n"
         << "D: "<< D<< "\n"
         << "Liczba kolumn: "<< IloscKolumn<< "\n"
         << "h: "<< h<< "\n"
         << "dt: "<< dt<< "\n"
         << "Liczba wierszy = "<< IloscWierszy<< "\n\n"
         << "________________________________________________________________________________\n";
}

/*metoda klasyczna bezposrednia*/
void KMetodaBezposrednia()
{
    double Temp_t = 0.0;

    for(int k = 1; k < IloscWierszy; k++)
    {
        for(int i = 1; i < IloscKolumn-1; i++)
        {
            Tablica[k][i]=Tablica[k-1][i]+Lambda*(Tablica[k-1][i-1]-(2*Tablica[k-1][i])+Tablica[k-1][i+1]);
        }
    }
    // Algorytm zoptymalizowany do podpunktu (1)

    /*for(int k = 0; k < IloscWierszy; k++)
    {
       double Temp_x = 0;
       double MaxBlad = 0.0;
       if(k == IloscWierszy-1)
       {
           for(int in = 0; in < IloscKolumn; in++)
           {
               if(fabs(Tablica[k][in] - RownanieAnalityczne(Temp_x,Temp_t)) > MaxBlad)
               {
                   MaxBlad = fabs(Tablica[k][in] - RownanieAnalityczne(Temp_x,Temp_t));
               }
           Temp_x += h;
           }
           Pomiary << "MaxBlad: \t" <<MaxBlad << endl;
       }
       Temp_t += dt;
   }*/


    // Algorytm zoptymalizowany do podpunktu (3)

    for(int k = 0; k < IloscWierszy; k++)
    {
        double Temp_x = 0;
        double MaxBlad = 0.0;
        for(int in = 0; in < IloscKolumn; in++)
        {
            if(fabs(Tablica[k][in] - RownanieAnalityczne(Temp_x,Temp_t)) > MaxBlad)
            {
                MaxBlad = fabs(Tablica[k][in] - RownanieAnalityczne(Temp_x,Temp_t));
            }
            Temp_x += h;
        }
        MaxBledy << Temp_t << "\t" << MaxBlad << endl;
        if(k == IloscWierszy-1)
        {
            Pomiary << "MaxBlad: \t" <<MaxBlad << endl;
        }
        Temp_t += dt;
    }
}




/*metoda dyskretyzacji Laasonena*/
void Laasonen()
{
    double *Ni = new double[IloscKolumn];
    double *L = new double[IloscKolumn];
    double *D=new double[IloscKolumn];
    double *U=new double[IloscKolumn];
    double *WektorB=new double[IloscKolumn];
    double *WektorX=new double[IloscKolumn];

    /*ustalam poczatkowe t*/
    double Temp_t = 0.0;

    /*Tworze odpowiednie macierze L,D,U*/
    for(int i= 0; i< IloscKolumn; i++)
    {
        L[i]= Lambda;
        D[i]= -(1.0+ (2*Lambda));
        U[i]= Lambda;
    }

    /*Uzupelniam macierze L,D,U o zalozenia*/
    L[IloscKolumn-1]= 0.0;
    U[0]= 0.0;
    D[0]= 1.0;
    D[IloscKolumn-1]= 1.0;

    /*wypelniam wektor Ni, ktory jest pomocny przy rozwiazywaniu ukladu metoda Thomasa*/
    Ni = TworzNi(L,D,U);

    /*przechodze do dyskretyzacji metoda Laasonena*/
    for(int k= 0; k< IloscWierszy; k++)
    {
        /*ustalam brzegowe x*/
        double Temp_x=0.0;

        /*wpisuje ustalone warunki brzegowe*/
        WektorB[0]= 1.0;
        WektorB[IloscKolumn-1]= 0.0;

        /*w tym przypadku juz nie wystepuje mnozenie macierzy przez wektor*/
        for(int i= 1; i< IloscKolumn-1; i++)
        {
            WektorB[i] = -(Tablica[k][i]);
        }

        /*rozwiazujemy nasz uklad metoda Thomasa*/
        Thomas(L,Ni,U,WektorX,WektorB);

        /*Uzupelniam nasza macierz wynikowa o wyliczony wektor X*/
        if(k== IloscWierszy- 1);
        else
        {
            for(int m= 0; m< IloscKolumn; m++)
                Tablica[k+ 1][m]= WektorX[m];
        }

        /*ustalam blad maksymalny*/
        double MaxBlad= 0.0;
        for(int in= 0; in< IloscKolumn; in++)
        {
            if(fabs(WektorX[in] - RownanieAnalityczne(Temp_x,Temp_t)) > MaxBlad)
            {
                MaxBlad = fabs(WektorX[in] - RownanieAnalityczne(Temp_x,Temp_t));
            }
            Temp_x += h;
        }

        /*zapisuje blad maksymalny do pliku*/
        MaxBledy<< Temp_t<< "\t"<< MaxBlad<< "\n";

        /*Tutaj pomocniczo zapisuje sobie czas wykonania i blad dla t_max do oddzielnego pliku*/
        if(k== IloscWierszy- 1)
        {
            Pomiary<< "MaxBlad: \t"<< MaxBlad<< "\n";
        }
        /*zwiekszam t o przyrost czasowy dt*/
        Temp_t += dt;

    }
}

/*metoda sluzaca do zapisywania rezultatow obliczania analitycznego do pliku*/
void DrukujAnalityczna(int IloscWierszy, int IloscKolumn, double x, double h)
{
    Analitycznie.open("TablicaAnalityczna.txt");
    /*tworze sobie tablice*/
    double **Tablica = new double *[IloscWierszy];
    for(int i= 0; i< IloscWierszy ; i++)
    {
        Tablica[i] = new double [IloscKolumn];
    }

    /*Wpisuje sobie najpierw arugenty do pliku*/
    Analitycznie << "\t";
    for(int k = 0; k < IloscKolumn ; k++)
    {
        Analitycznie<< x<< "\t";
        x+= h;
    }
    Analitycznie<< "\n";


    /*teraz obliczam wyniki dla rownania analitycznego przy zmiennym czasie T i zapisuje do pliku*/

    double t= 0.0;
    for(int i= 0; i< IloscWierszy; i++)
    {
        double x= -a;
        Analitycznie<< t<< "\t";
        for(int k= 0; k< IloscKolumn; k++)
        {
            Analitycznie<< RownanieAnalityczne(x,t)<< "\t";
            x+= h;

        }
        Analitycznie<< "\n";
        t+= dt;
    }
}

/*metoda pomocnicza sluzaca do obliczania czasu (w sekundach) wykonywania sie programu*/
double Czas()
{
    struct timeb czas;
    double sekundy;
    ftime(&czas);

    sekundy= (double) czas.time;
    sekundy += (double) czas.millitm / 1000.0;
    return sekundy;
}

int main()
{
    /*otwieram sobie pliki z pomiarami i z bledami*/
    Pomiary.open("Pomiary.txt");
    MaxBledy.open("BledyMax.txt");
    /*uzupelniam pierwszy wiersz pliku z bledami*/
    MaxBledy<< "t"<< "\t"<< "MaxBlad"<< "\n";
    /*ilosc kolumn macierzy*/
    int N=500;
    /*zmienne pomocnicze przy liczeniu czasu*/
    double pocz,kon;

    /*inicjuje sobie dane dla naszego przypadku*/
    inicjujDane(N);

    /*wypisuje sobie dane na ekran*/
    WypiszDane();

    /*tworze wstepna tablice wynikowa*/
    Tablica=TworzMacierz(IloscWierszy,IloscKolumn);
    /*Uzupelniam tablice o warunki poczatkowe i brzegowe*/
    Warunki();

    /*uruchamiam licznik czasu*/
    pocz= Czas();

    /*wykonuje metode Laasonena lub Cranka-Nicholsona*/
    KMetodaBezposrednia();
    //Laasonen();

    /*teraz sprawdzam ile czasu uplynelo*/
    kon= Czas();

    /*zapisuje sobie wyniki analityczne do pliku*/
    DrukujAnalityczna(IloscWierszy,IloscKolumn,-a, h);

    /*Na koniec zapisuje sobie do pliku pomocniczego czas wykonywania, krok H i iloscKolumn N*/
    Pomiary<< "Czas[sek]: \t"<< (kon- pocz)<< "\n";
    Pomiary<< "H: \t"<< h<< "\n";
    Pomiary<< "N: \t"<< IloscKolumn<< "\n";
    Pomiary<< "\n";

    /*zapisuje sobie obliczona macierz do pliku*/
    DrukujMacierz(Tablica, IloscWierszy, IloscKolumn, -a, h, t0, dt);
    cout<< "Wszystkie wyniki zostaly pomyslnie zapisane do plikow...\n";
    getchar();
    Pomiary.close();
    MaxBledy.close();
    return 0;
}


