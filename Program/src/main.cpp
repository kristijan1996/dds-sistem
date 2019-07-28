/**
 * @file main.c
 * @brief Implementacija DDS sistema
 *
 * U ovom programu implementiran je DDS sistem za generisanje
 * sinusnog signala proizvoljne ucestanosti iz opsega [0, 40] MHz.
 * Sve softverske instance koje imitiraju hardverske komponente su 
 * ogranicene fiksnom aritmetikom koriscenjem "ac_int.h" i "ac_fixed.h"
 * biblioteka Mentor Graphics Corporation-a. 
 * Korisnik zadaje zeljenu ucestanost, pocetnu fazu i velicinu ditera.
 * Taktovanje sistema se simulira velikom FOR petljom prolaskom kroz
 * koju se fazni akumulator inkrementira, CORDIC algoritam generise
 * sinusni i kosinusni signal koji se potom filtriraju FIR filtrom i 
 * na koje se dodaje diter koji razbija spurove u spektru. 
 * Neki medjurezultati se loguju u fajl koji MATLAB koristi za iscrtavanje
 * signala i spektralnu analizu dobijenih rezultata.
 * 
 * @date Jun 2019
 * @author Kristijan Mitrovic (mk150214d@etf.bg.ac.rs)
 * @author Dragan Bozinovic (bd150211d@etf.bg.ac.rs)
 *
 */

#include <iostream>
#include <fstream>
#include <random>
#include "ac_int.h"
#include "ac_fixed.h"
using namespace std;

/* Definicije makroa */
#define PI 3.141592
#define FIR_ORDER 8

int main()
{
    /* Korisnicki interfejs: frekvencija u Hz u opsegu [0, 40MHz],
    pocetna faza u opsegu [0, 2*PI] i ditherAmp koji diktira velicinu
    ditera koji se dodaje na fazni akumulator */
    float f0 = 1000000;  // Preporucene frekvencije su oblika fs/2^x
    float phi0 = 0;
    int ditherAmp = 3;

    /* W je broj bita kontrolne reci f0 i faznog akumulatora.
    N=14 zbog DAC-a, a M=14 je najveca smislena vrednost za ulaz u generator odbiraka*/
    const int W = 40, N = 14, M = N;
    /* Ucestanost odabiranja */
    float fs = 100000000.0;

    /* Alocira se fazni akumulator */
    ac_int<W, false> phaseAcc = 0; 
    /* dteta je frekvencija prevedena u inkrement akumulatora */
    ac_int<W, false> dteta = f0/fs*pow(2,W);
    /* phi je pocetna faza izrazena u inkrementu akumulatora */
    ac_int<W, false> phi = phi0 / (2*PI) * pow(2,W);

    /* Akumulator se inicijalizuje pocetnom fazom */
    phaseAcc = phi;

    /* Opseg [0, 2^(M-2)] mapiran je u opseg [0, PI/2] tako sto je u
    phaseAccOut izdvojeno M-2=12 nizih bita faznog akumulatora, */
    ac_int<M-2, false> phaseAccOut;
    /* a dva najvisa bita se koriste za odredjivanje kvadranta u kom se ugao nalazi  */
    int muxSel; 

    /* Promenljive za smestanje izlaznih odbiraka CORDIC-a */
    ac_fixed<N, 2, true> sampleSin, sampleCos;
    /* Nizovi za smestanje poslednjih FIR_ORDER+1 izlaznih odbiraka CORDIC-a koji
    idu na filtriranje */
    ac_fixed<N, 2, true> cordicSin[FIR_ORDER+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    ac_fixed<N, 2, true> cordicCos[FIR_ORDER+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    /* Koeficijenti FIR filtra za kompenzaciju sinx/x */ 
    ac_fixed<N, 2, true> filter[FIR_ORDER+1] = {0.00170522, -0.00583712,  0.01786389, -0.06833815, 0.81251125,
                                                 -0.06833815, 0.01786389, -0.00583712,  0.00170522};
    /* Promenljive za smestanje izlaza filtra, tj. ulaza u DAC */
    ac_fixed<N, 2, true> DACinputSin, DACinputCos;

    /* Popunjavanje tabele potrebne za CORDIC pre pocetka rada sistema. Posto je arctg u opsegu
    (-PI/2, PI/2), dovoljna su dva bita za ceo deo, od kojih je jedan za znak. */    
    ac_fixed<M, 2, true> atanAngle[M];
    for (int i = 0; i < M; i++)
    {
        atanAngle[i] = atan(pow(2,-i));
    }

    /* Instanciranje generatora nasumicnih brojeva normalne raspodele N ~ (0,1) */
    default_random_engine generator;
    normal_distribution<double> distribution(0,1);
    /* dither = ditherAmp*distribution(generator)*pow(2,W-M) je vrednost koja se 
    dodaje na fazni akumulator */
    double dither;

    /* Otvara se fajl za logovanje rezultata */
    ofstream file;
    file.open("DDSlog.csv");
    /* Ukoliko je neuspesno otvaranje fajla, iskociti iz programa */
    if (!file.is_open())
    {
        cout << "Neuspesno otvaranje fajla za logovanje!" << endl;
        return -1;
    }

    /* Zapis frekvencije sinusoide u fajl zarad lakseg iscrtavanja */
    file << f0 << "\n";

    /* Velika FOR petlja koja simulira taktovanje sistema. U fajl ce biti
    izlogovano 2^16*dt [sec] rada sistema, gde je dt period izmedju dva takta */
    for (int idx = pow(2,16); idx > 0; idx--)
    {
        /* 12 bita se uzima za mapiranje opsega 0-pi/2 u 0-2^12 */
        phaseAccOut = phaseAcc.slc<M-2>(W-M);
        /* a preostala dva za odredjivanje kvadratna u kom se ugao nalazi */
        muxSel = phaseAcc.slc<2>(W-2).to_int();

        /* Inicijalizacija CORDIC algoritma */
        ac_fixed<N, 2, true> x = 0.60725294, y = 0, x_new = 0, y_new = 0;
        /* Brojna vrednost iz faznog akumulatora se skalira u ugao u opsegu [0, PI/2] */
        ac_fixed<N, 2, true> z = phaseAccOut.to_double()/(pow(2,M-2)-1) * PI/2; 
        int sigma = (z > 0 ? 1 : -1);
        /* Broj iteracija CORDIC-a jednak je duzini binarne reci */
        for (int i = 0; i < M; i++)
        {
            x_new = x - sigma*(y >> i);
            y_new = y + sigma*(x >> i);
            z = z - sigma*atanAngle[i];
            sigma = (z > 0 ? 1 : -1);  

            x = x_new; y = y_new;       
        }

        /* Zapis sadrzaja faznog akumulatora u fajl */
        file << phaseAcc << ", " << phaseAccOut << ", ";

        /* muxSel nosi informaciju o kvadrantu u kom se ugao nalazi, na osnovu cega
        se odredjuje konacan izlaz CORDIC-a */
        switch(muxSel)
        {
            case 0:
                sampleSin = y; 
                sampleCos = x;
                /* Ispis u fajl */
                file << y << ", " << x <<  ", ";
                break;
            case 1:
                sampleSin = x; 
                sampleCos = -y;
                /* Ispis u fajl */
                file << x << ", " << -y << ", ";
                break;
            case 2:
                sampleSin = -y; 
                sampleCos = -x;
                /* Ispis u fajl */
                file << -y << ", " << -x << ", ";
                break;
            case 3:
                sampleSin = -x; 
                sampleCos = y;
                /* Ispis u fajl */
                file << -x << ", " << y << ", ";
                break;
        }

        /* Pamcenje poslednjih FIR_ORDER+1 izlaza CORDIC-a koji potom idu na filtriranje */
        for (int i=0; i<FIR_ORDER; i++)
        {
            cordicSin[i] = cordicSin[i+1];
            cordicCos[i] = cordicCos[i+1];
        }
        cordicSin[FIR_ORDER] = sampleSin;
        cordicCos[FIR_ORDER] = sampleCos;

        /* Filtriranje signala */
        DACinputCos = 0;
        DACinputSin = 0;
        for (int i=0; i<FIR_ORDER+1; i++)
        {
            DACinputCos += filter[i]*cordicCos[FIR_ORDER-i];
            DACinputSin += filter[i]*cordicSin[FIR_ORDER-i];  
        }

        /* Ispis u fajl */
        file << DACinputSin << ", " << DACinputCos << "\n";

        /* Dodavanje ditera za poptiskivanje spurova */
        dither = ditherAmp*distribution(generator)*pow(2,W-M);

        /* Na svaki takt se vrednost faznog akumulatora
        inkrementira za diterovano dteta */
        phaseAcc += (ac_int<W, false>)(dteta.to_double() + dither);
    }

    /* Zatvaranje fajla za logovanje */
    file.close();
    
    return 0;
}