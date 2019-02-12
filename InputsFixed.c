/* Projekt 2 - Iteracne Vypocty
 * Meno: Peter Vinarcik
 * Login: <xvinar00>
 * Email: <xvinar00@stud.fit.vutbr.cz>
 * Datum: 13.11.2018
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

double cfrac_log(double x, unsigned int n);                                         //deklaracia zretazeny zlomok_logartitmus
double taylor_log(double x, unsigned int n);                                        //deklaracia taylorov_polynom
double log_kontrola(double x);                                                      //deklaracia kontrola vypoctov <math.h>
double taylor_pow(double x, double y, unsigned int n);                              //deklaracia taylorov exp funkcia
double taylorcf_pow(double x, double y, unsigned int n);                            //deklaracia zretaz_zlomok exp funkcia
double pow_kontrola(double x, double y);                                            //deklaracia kontrola vypoctov <math.h>
int arg_kontrola();                                                                 //deklaracia kontrola argv[1]
double mylog(double x);                                                             //deklaracia pre zistovanie efektivnejsieho logaritmu
double mypow(double x, double y);                                                   //deklaracia pre mypow -> vlastne mocnenie

int main(int argc, char *argv[])
{
    if (argc == 4){                                                                 //4 argumenty pre --log
        if (strcmp(argv[1],"--log") != 0) {                                         //argv[1] sa musi rvonat --log
            arg_kontrola();                                                         //stderr -> chybove hlasenie
        }else{
            char *znakN, *znakX;                                                    //premenne na kontrolu vstupu
            double x = strtod(argv[2], &znakX);                                     //nacitanie x, strtod (string to double)
            unsigned long int n = strtoul(argv[3], &znakN, 10);                     //strtoul poc. iteracii n to u_l_int
            if (*znakN != '\0' || *znakX != '\0' || n == 0 || argv[3][0] == '-'){
                arg_kontrola();
                return EXIT_FAILURE;
            }
            log_kontrola(x);                                                        //kontrola logaritmu z <math.h>
            printf(" cfrac_log(%.5g) = %.12g\n",x, cfrac_log(x,n));                 //vypise cflog na 12 des. miest
            printf("taylor_log(%.5g) = %.12g\n", x, taylor_log(x,n));               //vypis taylor_log pre dany interval resp. podmienku
            printf("     mylog(%.5g) = %.7e\n",x, mylog(x));                        //ak bola presnejsia taylog vypise <--
            return 0;
        }
    }
    else if (argc == 5){                                                            //5 argumentov pre --pow
        if (strcmp(argv[1],"--pow") != 0) {                                         //argv[1] sa musi rovnat --pow
            arg_kontrola();                                                         //stderr -> chybove hlasenie
        }else{
            char *znakN, *znakX, *znakY;                                            //premenne na kontrolu vstupu
            double x = strtod(argv[2], &znakY);                                     //nacitanie x, strtod (string to double)
            double y = strtod(argv[3], &znakX);                                     //nacitanie x, strtod (string to double)
            unsigned long int n = strtoul(argv[4], &znakN, 10);                     //strtoul poc. iteracii n to u_l_int
            if (*znakN != '\0' || *znakX != '\0' || *znakY != '\0' || n == 0 || argv[4][0] == '-'){
                arg_kontrola();                                                     //podmienka ak vstup neobsahuje iba cislo -> EXIT_FAILURE
                return EXIT_FAILURE;
            }
            pow_kontrola(y, x);                                                     //kontrola umocenenia z <math.h>
            printf("  taylor_pow(%.5g,%.5g) = %.12g\n", x, y, taylor_pow(x,y,n));   //vypis exp. funkcie pomocou talyorovho polynomu
            printf("taylorcf_pow(%.5g,%.5g) = %.12g\n", x, y, taylorcf_pow(x,y,n)); //vypis exp. funkcie pomocou talyorovho polynomu
            printf("       mypow(%.5g,%.5g) = %.7e\n",x, y, mypow(x,y));
            return 0;
        }
    } else{
        arg_kontrola();                                                             //EXIT_FAILURE v pripade ze pocet arugmentov nieje 4/5
    }
}

double cfrac_log(double x, unsigned int n)                                          //funkcia na vypocet cfrac_log
{
    if (x == 0){
        return -INFINITY;
    }
    if(isinf(x)){
        return INFINITY;
    }
    double cflog = 0.0;                                                             //vysledok -> ulozi sa do tejto prem.
    double c = 0;                                                                   //c = 0 --> start na dne [0]
    double z = (x - 1) / (x + 1);                                                   //vyjadrenie z zo vzorca
    for (int i = n; i >= 2; i--){                                                   //cyklus pre iteracie
        c = ((n - 1) * (n - 1) * z * z) / ((2*n - 1) - c);                          //vzorec na vypocet c, mimo 1-c
        n--;                                                                        //znizenie n, pretoze ideme z dna
    }
    if (n == 1){                                                                    //ked n=1 --> posledna iteracie (delenie)
        cflog = 2*z / (1-c);                                                        //vypocet cflog
    }
    return cflog;                                                                   //vracia vysleddok cflog
}

double taylor_log(double x, unsigned int n)                                         //funkcia na vypocet tay_log
{
    if (x == 0){
        return -INFINITY;
    }
    if(isinf(x)){
        return INFINITY;
    }
    double taylog = 0.0;
    double nove_y = 1.0;
    if (x<1){                                                                       //pre interval x E(0,1)
        double y = 1 - x;                                                           //substitucia pretoze vzorec je pre (1-y)
        for (unsigned int i = 1; i <= n; i++) {                                     //opakovanie cyklu podla pocut iteracii
            nove_y *= y;                                                            //n_y = n_y * y (zacina na 1.0 kvoli prvej iteracii)
            taylog -= nove_y / i;                                                   // - -> podla vzorca -y/1 - y^2/2
        }
    }else{                                                                          //pre interval x E (1 az nekonecno)
        for (unsigned int i = 1; i <= n; i++){                                      //opakovanie podla poctu iteracii
            nove_y *= ((x - 1) / x);                                                //vzorec na vypocet ((x-1 / x)^n zacina na 1 kvoli prvej iteracii
            taylog += (nove_y / i);                                                 //kvoli EPS --> nove_y / i   potom *i++*
        }
    }
    return taylog;
}

double log_kontrola(double x)                                                       //porovnanie s matematickou kniznicou
{
    printf("       log(%.5g) = %.12g\n",x, (log(x)));                               //logaritmus z <math.h> + -lm
    return log(x);
}

double taylor_pow(double x, double y, unsigned int n)                               //funkcia na vypocet exponenc pomocou taylor_polynom
{
    double fakt = 1.0;
    double tl_log;                                                                  //premenna na priradenie logaritmu z taylorovho polynomu
    tl_log = taylor_log(x, n);                                                      //priradenie log
    double vys = 1.0, new_ln = 1.0, new_y = 1.0;
    double celk_taylor = 1.0;
    for (unsigned int i = 1; i < n; i++) {                                          //cyklus opakovanie podla poctu iteracii (vzdy n-1 iteracii lebo prva je 1)
        fakt *= i;                                                                  //vzorec na faktorial
        new_y *= y;                                                                 //x, x^2, x^3...
        new_ln *= tl_log;                                                           //ln a,ln^2a (lna * lna)...
        vys += (new_y * new_ln)/fakt;                                               //vysledok podla vzorca
        if (isinf(vys) || isnan(vys))                                               //zisti ci cislo v urcitej iteracii nevyhodnoti ako inf, vrati predchadzajucu iteraciu
            break;
        celk_taylor = vys;
    }
    return celk_taylor;
}

double taylorcf_pow(double x, double y, unsigned int n)                             //funkcia na vypocet exponenc pomocou zretazeneho zlomku
{
    double fakt = 1.0;
    double cf_log;                                                                  //premenna na priradenie logaritmu z taylorovho polynomu
    cf_log = cfrac_log(x, n);                                                       //priradenie log
    double vys2 = 1.0, new_ln = 1.0, new_y = 1.0;
    double celk_cf = 1.0;
    for (unsigned int i = 1; i < n; i++) {                                          //cyklus opakovanie podla poctu iteracii (vzdy n-1 iteracii lebo prva je 1)
        fakt *= i;                                                                  //vzorec na faktorial
        new_y *= y;                                                                 //x, x^2, x^3...
        new_ln *= cf_log;                                                           //ln a,ln^2a (lna * lna)...
        vys2 += (new_y * new_ln)/fakt;                                              //vysledok podla vzorca
        if (isinf(vys2) || isnan(vys2))                                             //zisti ci cislo v urcitej iteracii nevyhodnoti ako inf, vrati predchadzajucu iteraciu
            break;
        celk_cf = vys2;
    }
    return celk_cf;
}

double pow_kontrola(double x, double y)                                             //porovnanie s matematickou kniznicou
{
    printf("         pow(%.5g,%.5g) = %.12g\n",y, x, (pow(y,x)));
    return pow(y,x);                                                                //vypis + umocnenie z <math.h> + -lm <pre kontrolu>
}

int arg_kontrola()                                                                  //stderr --> kontrola argv[1]
{
    fprintf(stderr,"Zle zadane argumenty.\n");
    return EXIT_FAILURE;
}

double mylog(double x)
{                                                                                   //VZOREC PRE VYPOCET TAYLOROV POLYNOM!!
    double my_taylog = 0.0, my_cflog = 0.0;
    double nove_y = 1.0;
    double rozdiel_taylog = 1.0, rozdiel_cflog = 1.0;                               //novy vysledok - stary vysledok = rozdiel
    double taylog_old = 0.0, cflog_old = 0.0;                                       //drzanie predosleho vysledku kvoli porovnaniu presnosti
    double eps = 0.00000001;                                                        //presnost na %.7e
    int poc_taylog = 0, poc_cflog = 0;                                              //premenne na pocitanie itearcii
    unsigned int i = 1, n = 1;
    while (rozdiel_taylog > eps){
        if (x < 1) {                                                                //pre interval x E(0,1)
            double y = 1 - x;                                                       //substitucia pretoze vzorec je pre (1-y)
            nove_y *= y;                                                            //n_y = n_y * y (zacina na 1.0 kvoli prvej iteracii)
            my_taylog -= nove_y / i;                                                // - -> podla vzorca -y/1 - y^2/2
            rozdiel_taylog = (fabs(my_taylog) - fabs(taylog_old));
            taylog_old = my_taylog;                                                 //taylog sa na dalsi cyklus zmeni, stary bude ulozeny v old
            poc_taylog++;                                                           //pocitanie poctu iteracii
            i++;
        } else {                                                                    //pre interval x E (1 az nekonecno)
            nove_y *= ((x - 1) / x);                                                //vzorec na vypocet ((x-1 / x)^n zacina na 1 kvoli prvej iteracii
            my_taylog += (nove_y / i);                                              //kvoli EPS --> nove_y / i   potom *i++*
            rozdiel_taylog = (fabs(my_taylog) - fabs(taylog_old));
            taylog_old = my_taylog;                                                 //taylog sa na dalsi cyklus zmeni, stary bude ulozeny v old
            poc_taylog++;                                                           //pocitanie poctu iteracii
            i++;
        }
    }                                                                               //VZOREC PRE VYPOCET ZRETAZENEHO ZLOMKU!!!
    while(rozdiel_cflog > eps){                                                     //KYM JE ROZDIEL VACSI AKO EPS BUDE CYKLUS AKTIVNY
        unsigned int pom;                                                           //pomocna ktora bude drzat nasledujucu hodnotu N
        pom = n+1;                                                                  //n+1 pretoze najprv vykona iteraciu "poslednu" a potom dalsiu, porovna 1 s 2, 2 s 3 atd. atd.
        double c = 0;                                                               //c = 0 --> start na dne [0]
        double z = (x - 1) / (x + 1);                                               //vyjadrenie z zo vzorca
        while (n != 1) {                                                            //ked sa n == 1 dokonci vypocet citatel / menovatel
            c = ((n - 1) * (n - 1) * z * z) / ((2 * n - 1) - c);
            n--;                                                                    //vzorec na vypocet c, mimo 1-c
        }                                                                           //n-- aby sme sa dostali na podmienku n == 1
        my_cflog = 2 * z / (1 - c);                                                 //vypocet my_cflog
        rozdiel_cflog = (fabs(my_cflog) - fabs(cflog_old));                         //odcitanie stareho a noveho, ak eps = 0.00000001 cyklus konci
        cflog_old = my_cflog;                                                       //mycflog sa zmeni na dalsej iteracii,stary osatne v old
        poc_cflog++;                                                                //pocitanie poctu iteracii
        n = pom;                                                                    //predtym bolo pred 1 delenie, teraz pre 2, pootm pre 3, 4, 5 atd. atd.
    }

    if (poc_taylog < poc_cflog){                                                    //zisti, ktora funkcia mala mensi pocet iteracii kedze cyklus skoncil
        return my_taylog;
    }else{
        return my_cflog;
    }
}

double mypow(double x, double y)                                                    //funkcia pre umocnovanie pomocou efektivnejsieho logaritmu
{
    double fakt = 1.0;
    double cf_log;                                                                  //premenna na priradenie logaritmu z mylog
    cf_log = mylog(x);                                                              //priradenie log
    double my_pow = 1.0, new_ln = 1.0, new_y = 1.0;                                 // f(x) ^ 2 + my_pow pre prvu iteraciu 1+...
    double rozdiel_pow = 1.0;                                                       //podmienka pre while (zacina na 1.0 aby sa while spustil, neskor sa prepise
    double old_pow = 0.0;                                                           //premenna pre stary vysledok
    unsigned int i = 1;
    double eps = 0.00000001;                                                        //presnost eps
    while (rozdiel_pow > eps){
        fakt *= i;                                                                  //vzorec na faktorial
        new_y *= y;                                                                 //x, x^2, x^3...
        new_ln *= cf_log;                                                           //ln a,ln^2a (lna * lna)...
        my_pow += (new_y * new_ln)/fakt;                                            //vysledok podla vzorca
        rozdiel_pow = (fabs(my_pow) - fabs(old_pow));                               //zisti ci stale plati ze rozdiel > eps
        old_pow = my_pow;                                                           //priradi predosli vysledok do old_pow na porovnanie v dalsej iteracii
        i++;                                                                        //zvacsovanie i pre faktorial
    }
    return my_pow;
}