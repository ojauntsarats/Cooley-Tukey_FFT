#include <iostream>
#include <vector>
#include <math.h>
#include <complex>

double PI = 3.14159265358979;

using namespace std;

bool check_timePoints(int _size){
    //Compruebo si el número de puntos es una potencia de dos
    return ((_size !=0) && ((_size & (_size-1)) == 0));
}

void Radit2_DIT(std::complex<double> X[],  int N){

    if(N==1){
        return;
    }
    else{

        //divide
        std::complex<double> E[N/2];
        std::complex<double> O[N/2];
        for(int i=0; i<N/2; i++){
            E[i] = X[2*i];
            O[i] = X[2*i+1];
        }

        Radit2_DIT(E, N/2);
        Radit2_DIT(O, N/2);
        for(int k=0; k<N/2; k++){
            X[k] = E[k] + O[k]*std::polar(1.0, -2*PI*k/N);
            X[k + N/2] = E[k] - O[k]*std::polar(1.0, -2*PI*k/N);
        }
    }
}

std::vector<int> factorize(int _number){
    std::vector<int> a;
    int i=2;
    int j=0;

    while(_number>1){
        if(_number%i==0){
            _number=_number/i;
            a.push_back(i);
            j++;
            i=2;
        }
        else
            i++;
    }
    return a;
}

void DFT_complex1D(complex<double> X[],int N){

    std::cout << "Performing 1D Complex DFT..." << std::endl;

        unsigned int nTimeSteps = N;
        std::vector<std::complex<double>> x;
        x.resize(nTimeSteps);

        for(unsigned int iTime=0; iTime<nTimeSteps; iTime++){
            x[iTime] = X[iTime];
        }

        std::vector<std::complex<double>> _X;
        _X.resize(nTimeSteps);


        for(unsigned short p=0; p<nTimeSteps; p++){
            for(unsigned short k=0; k<nTimeSteps; k++){
                std::complex<double> expI = -2*1i*PI*p*k/nTimeSteps;
                _X[p] += x[k]*std::exp(expI);
            }
        }

        for(int i=0; i<N; i++){
            X[i] = _X[i];
        }

}

void Cooley_Tukey_FFT(complex<double> X[], int N){


    if(check_timePoints(N)){
        Radit2_DIT(X, N);
    }
    else if(factorize(N).size()==1){
        DFT_complex1D(X, factorize(N)[0]);
    }
    else{

        int _n1 = factorize(N)[factorize(N).size()-1];
        int _n2 = N/_n1;

        complex<double> matrix[_n2][_n1];
        for(int i=0; i<_n2; i++){
            for(int j=0; j<_n1; j++){
                matrix[i][j] = X[_n1*i+j];
            }
        }

        for(int j=0; j<_n1; j++){
            complex<double> dft_in[_n2];
            for(int i=0; i<_n2; i++){dft_in[i] = matrix[i][j];}
            if(check_timePoints(_n2)){Radit2_DIT(dft_in, _n2);}
            else{Cooley_Tukey_FFT(dft_in,_n2);}
            for(int i=0; i<_n2; i++){ matrix[i][j] = dft_in[i];}
        }

    complex<double> factor[_n2][_n1];
    for(int b=0; b<_n2; b++){
        for(int j=0; j<_n1; j++){
            factor[b][j] = polar(1.0,-2*PI*b*j/N);
        }
    }

    for(int b=0; b<_n2; b++){
        for(int j=0; j<_n1; j++){
            matrix[b][j] = matrix[b][j]*factor[b][j];
        }
    }

    for(int b=0; b<_n2; b++){
        complex<double> dft_out[_n1];
        for(int a=0; a<_n1; a++){dft_out[a] = matrix[b][a];}
        if(check_timePoints(_n1)){Radit2_DIT(dft_out,_n1);}
        else{Cooley_Tukey_FFT(dft_out,_n1);}
        for(int a=0; a<_n1; a++){matrix[b][a] = dft_out[a];}
    }

    for(int a=0; a<_n1; a++){
        for(int b=0; b<_n2; b++){
            X[_n2*a+b] = matrix[b][a];
        }
    }
    }
}

int main()
{
    std::vector<double> x = {0.035063,0.0349655,0.0348717}; //,0.0347853,0.0347095,0.0346473,0.0346011,0.0345726}; //0.034563,0.0345726,0.034601,0.0346472,0.0347094,0.0347852,0.0348716,0.0349654,0.0350629,0.0351604,0.0352542,0.0353406,0.0354164,0.0354786,0.0355248,0.0355533, 0.0355629,0.0355533,0.0355249,0.0354787,0.0354166,0.0353409}; //,0.0352545, 0.0351607,0.0350632, 0.0349657 ,0.0348719};

    complex<double> x_t[x.size()];

    for(unsigned int i=0; i<x.size(); i++){
        x_t[i].real(x[i]);
        x_t[i].imag(0.0);
    }

    Cooley_Tukey_FFT(x_t, x.size());

    for(unsigned int i=0; i<x.size(); i++){
        cout << sqrt(x_t[i].real()*x_t[i].real() + x_t[i].imag()*x_t[i].imag()) << endl;
    }


    return 0;
}
