#include "macierz.h"
#include "calerf.h"
#include <cstdio>
#include <fstream>
#include <sys/timeb.h>
#include <iostream>
#include <math.h>
#include <cstring>
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
    Lambda= 1;
    /*ustalamy t poczatkowe*/
    t0= 0.0;
    /*ustalamy t koncowe*/
    tMax=2.0;
    /*ustalamy parametry b i D*/
   // b= 0.1;
    D= 1.0;

    /*obliczamy a, h, t, dt oraz liczbe wierszy macierzy*/
    a= 6* sqrt( D* tMax );
    h= (a+ a)/ (IloscKolumn - 1);
    dt= (Lambda* h* h)/ D;
    d= (tMax/ dt)+ 1;
    IloscWierszy= (int)d;
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

void swap(double *W1,double *W2,int size){
    for (int j = 0; j < size; j++) {
        double tmp = W1[j];
        W1[j] = W2[j];
        W2[j] = tmp;
    }
}

void dekompozycja_LU(double **A, double *b, double *wyn)
{

    //Dekompozycja LU macierzy- metoda eliminacji Gaussa
    double x;
    int M=IloscKolumn;
    for(int k=0; k<M-1; k++)
    {
        for(int i=k+1; i<M; i++)
        {
            x = A[i][k]/A[k][k];
            A[i][k] = x;
            for(int j=k+1; j<M; j++)
            {
                A[i][j] = A[i][j] - (x*A[k][j]);
            }
        }
    }
    //Rozwi�zywanie uk�adu r�wna�
    double suma;
    double *z = new double[M];

    //podstawianie w prz�d
    for(int i=0; i<M; i++)
    {
        suma = 0;
        for(int j=0; j<=i-1; j++)
        {
            suma += A[i][j]*z[j];
        }

        z[i] = b[i]-suma;
    }

    //podstawianie w ty�
    for(int i=M-1; i>=0; i--)
    {
        suma = 0;
        for(int j=i+1; j<M; j++)
        {
            suma +=A[i][j]*wyn[j];
        }

        wyn[i] = (z[i]-suma)/A[i][i];
    }

    /*int size=IloscKolumn;
    double L[size][size];
    int permutation[size];


    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            L[j][i]=0;
        }
    }


    for(int i=0;i<size;i++){
        int indexMax = i;
        double MAX = A[i][i];

        //Wyznaczenie pivota
        for (int j = i; j < size; j++) {
            if (fabs(A[j][i]) > MAX) {
                MAX = fabs(A[j][i]);
                indexMax = j;
            }
        }
        permutation[i] = indexMax;

        //Jeśli znaleźlismy pivot to zamieniamy wiersz 1 z pivotem
        if (indexMax != i) {
            swap(A[indexMax], A[i], size);
            swap(L[indexMax], L[i], size);
        }
        //Tworzenie macierzy L i przeksztalcamy A
        for(int j=i+1;j<size;j++){
            double m = A[j][i] / A[i][i];
            L[j][i]=m;
            for(int k=i;k<size;k++) {
                A[j][k] -= (A[i][k] * L[j][i]);
            }
        }

        for(int j=0;j<size;j++){
            if(i==j) L[j][i]=1;
        }
    }

    for(int i=0;i<size;i++){
        if(i<size-1) {
            double tmp = b[i];
            b[i] = b[permutation[i]];
            b[permutation[i]] = tmp;
        }
    }

    for(int i=1;i<size;i++){
        for(int j=0;j<i;j++){
            b[i]-=L[i][j]*b[j];
        }
    }

    double y[size];
    for(int i=0;i<size;i++){
        y[i]=b[i];
        for(int j=0;j<i;j++){
            y[i]-=L[j][i]*y[j];
        }
    }

    for (int i = size-1; i >= 0; --i) {
        x[i]=y[i];
        for(int j=(i+1);j<size;j++){
            x[i]-=A[i][j]*x[j];
        }
        x[i]/=A[i][i];
    }*/
}

void DrukujPktAnal(double t, double x){
    string name= "t="+to_string(t)+"Analitycznie.csv";
    plik.open(name);
    plik << "\t";
    for(int k = 0; k < IloscKolumn ; k++)
    {
        plik << x << ";" << RownanieAnalityczne(x,t) << ";" << endl;
        x +=h;
    }plik << endl;

    plik.close();
}

void DrukujPktNum(double t, double x,double* Tablica){
    string name= "t="+to_string(t)+"Numerycznie.csv";
    plik.open(name);
    plik << "\t";
    for(int k = 0; k < IloscKolumn ; k++)
    {
        plik << x << ";" << Tablica[k] << ";" << endl;
        x +=h;
    }plik << endl;

    plik.close();
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
        MaxBledy << Temp_t << ";" << MaxBlad << endl;
        if(k == IloscWierszy-1)
        {
            Pomiary << "MaxBlad: ;" <<MaxBlad << endl;
        }

        if(k==IloscWierszy/3 || k==IloscWierszy-1){
            DrukujPktAnal(Temp_t,0.0);
            DrukujPktNum(Temp_t,0.0,Tablica[k]);
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
    double **A=new double*[IloscKolumn];
    for (int j = 0; j < IloscKolumn; ++j) {
        A[j]=new double[IloscKolumn];
        for (int i = 0; i < IloscKolumn; ++i) {
            A[j][i]=0;
        }
    }
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

    for (int l = 0; l < IloscKolumn; ++l) {
        try {
            A[l][l - 1] = L[l];
            A[l][l] = D[l];
            A[l][l + 1] = U[l];
        }
        catch(exception){
            continue;
        }
    }
    /*Uzupelniam macierze L,D,U o zalozenia*/
    A[IloscKolumn-1][IloscKolumn-1]= 0.0;
    A[0][1]= 0.0;
    A[0][0]= 1.0;
    A[IloscKolumn-1][IloscKolumn-1]= 1.0;

    double **tmp_A=new double*[IloscKolumn];
    for (int j = 0; j < IloscKolumn; ++j) {
        tmp_A[j] = new double[IloscKolumn];
    }

    /*przechodze do dyskretyzacji metoda Laasonena*/
    for(int k= 0; k< IloscWierszy; k++)
    {

        for (int j = 0; j < IloscKolumn; ++j) {
            for (int i = 0; i < IloscKolumn; ++i) {
                tmp_A[j][i]=A[j][i];
            }
        }

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
        //Thomas(L,Ni,U,WektorX,WektorB);
        dekompozycja_LU(tmp_A,WektorB,WektorX);

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
        MaxBledy<< Temp_t<< ";"<< MaxBlad<< ";\n";

        /*Tutaj pomocniczo zapisuje sobie czas wykonania i blad dla t_max do oddzielnego pliku*/
        if(k== IloscWierszy- 1)
        {
            Pomiary<< "MaxBlad: "<< MaxBlad<< ";\n";
        }
        /*zwiekszam t o przyrost czasowy dt*/
        cout << k << endl;

        if(k==IloscWierszy/3 || k==IloscWierszy-1){
            DrukujPktAnal(Temp_t,0.0);
            DrukujPktNum(Temp_t,0.0,Tablica[k]);
        }

        Temp_t += dt;
    }
}


/*metoda sluzaca do zapisywania rezultatow obliczania analitycznego do pliku*/
void DrukujAnalityczna(int IloscWierszy, int IloscKolumn, double x, double h)
{
    Analitycznie.open("TablicaAnalityczna.csv");
    /*tworze sobie tablice*/
    double **Tablica = new double *[IloscWierszy];
    for(int i= 0; i< IloscWierszy ; i++)
    {
        Tablica[i] = new double [IloscKolumn];
    }

    /*Wpisuje sobie najpierw arugenty do pliku*/
    Analitycznie << ";";
    for(int k = 0; k < IloscKolumn ; k++)
    {
        Analitycznie<< x<< ";";
        x+= h;
    }
    Analitycznie<< ";\n";


    /*teraz obliczam wyniki dla rownania analitycznego przy zmiennym czasie T i zapisuje do pliku*/

    double t= 0.0;
    for(int i= 0; i< IloscWierszy; i++)
    {
        double x= 0.0;
        Analitycznie<< t<< ";";
        for(int k= 0; k< IloscKolumn; k++)
        {
            Analitycznie<< RownanieAnalityczne(x,t)<< ";";
            x+= h;

        }
        Analitycznie<< ";\n";
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
    Pomiary.open("Pomiary.csv");
    MaxBledy.open("BledyMax.csv");
    /*uzupelniam pierwszy wiersz pliku z bledami*/
    MaxBledy<< "t"<< ";"<< "MaxBlad"<< ";\n";
    /*ilosc kolumn macierzy*/
    int N=100;
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

    /*wykonuje metode*/
    //KMetodaBezposrednia();
    Laasonen();

    /*teraz sprawdzam ile czasu uplynelo*/
    kon= Czas();

    /*zapisuje sobie wyniki analityczne do pliku*/
    DrukujAnalityczna(IloscWierszy,IloscKolumn,0, h);

    /*Na koniec zapisuje sobie do pliku pomocniczego czas wykonywania, krok H i iloscKolumn N*/
    Pomiary<< "Czas[sek]: \t"<< (kon- pocz)<< "\n";
    Pomiary<< "H: \t"<< h<< "\n";
    Pomiary<< "N: \t"<< IloscKolumn<< "\n";
    Pomiary<< "\n";

    /*zapisuje sobie obliczona macierz do pliku*/
    DrukujMacierz(Tablica, IloscWierszy, IloscKolumn, 0, h, t0, dt);
    cout<< "Wszystkie wyniki zostaly pomyslnie zapisane do plikow...\n";
    getchar();
    Pomiary.close();
    MaxBledy.close();
    return 0;
}


