// g++ ellipsoids.cpp -fopenmp -o ellipsoids

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 4
#define PI 3.14159265
#define INF 10.0
#define ppmin 0.5
#define ppmax 0.96
#define kappaoff 0.55

using namespace std;


float sqr(float x)
{
	return x*x;
}



int main(int args, char **argv)
{
	if(args != 2) {
		cout << "Usage: lens [delta]" << "\n";
		return 0;
	}

	const float delta = stof(argv[1]);
	cout << "\ndelta: " << delta << "\n";
	
	float worstalpha, worstsigma, worstpp, worstchi;
	float minvalue = INF;
	
	cout << "Max value of pp = " << ppmax << ", so max value of kappa = " << sqrt(2.0*ppmax*(sqrt(ppmax+2.0) - 1.0) - kappaoff) << ", ratio = " << (sqrt(2.0*ppmax*(sqrt(ppmax+2.0) - 1.0) - kappaoff))/(sqrt(0.6)/2.0) << "\n\n";
	
	const float alphamax = atan(sqrt(ppmax*(ppmax+1.0)));
	const float alphastep = alphamax / ceil(alphamax / (2.0*delta));
	const float sigmastep = PI / ceil(PI / (2.0*delta));;
	const float ppstep = (ppmax-ppmin) / ceil((ppmax-ppmin) / (2.0*delta));
	const float chistep = PI / (2.0 * ceil(PI / (4.0*delta)));
	
	double starttime = omp_get_wtime();
	
	int nthreads;
	omp_set_num_threads(NUM_THREADS);
	
	#pragma omp parallel
	{
		float alpha, sigma, pp, chi;
		float worstalphalocal, worstsigmalocal, worstpplocal, worstchilocal;
		float Cx, Cy, Dx, Dy, Sx, Sy, Zx, Zy, Xx, Xy, prl, prXx, prXy;
		float tanalpha, cosalphasigma, sinalphasigma;
		float value, minvaluelocal = INF;
		int id, nthr;
		
		id = omp_get_thread_num();
		nthr = omp_get_num_threads();
		if(id == 0)
		{
			nthreads = nthr;
			cout << "Number of threads = " << nthreads << "\n\n";
		}
		
		for(alpha = (id+0.5)*alphastep; alpha <= alphamax; alpha += nthreads*alphastep)
		{
			Cx = 2.0*sin(alpha);
			Cy = 1.0 - 2.0*cos(alpha);
			tanalpha = tan(alpha);
			
			for(sigma = sigmastep/2.0; sigma <= PI; sigma += sigmastep)
			{
				cosalphasigma = cos(alpha-sigma);
				sinalphasigma = sin(alpha-sigma);
				
				Sx = Cx - sinalphasigma;
				Sy = Cy + cosalphasigma;
				
				for(pp = ppmin + ppstep/2.0; pp <= ppmax; pp += ppstep)
				{
					for(chi = chistep/2.0; chi <= PI/2.0; chi += chistep)
					{
						Xx = sqrt(pp*(pp+1.0))*cos(chi);
						Xy = -pp*sin(chi);
						
						//if(Xx < tanalpha*(1.0-Xy)) break;
						
						prl = sqrt(sqr(Xx-Cx) + sqr(Xy-Cy));
						prXx = Cx + (Xx-Cx)/prl;
						prXy = Cy + (Xy-Cy)/prl;
						
						if((sqr(prXx)*pp + sqr(prXy)*(pp+1.0) >= pp*pp*(pp+1.0)) && (sqr(prXx-Sx) + sqr(prXy-Sy) <= 2.0*pp*(sqrt(pp+2.0)-1.0)-kappaoff))
						{
							value = pp*pp*(pp+1.0) - (sqr(cosalphasigma*(Xx-Sx) + sinalphasigma*(Xy-Sy))*pp + sqr(-sinalphasigma*(Xx-Sx) + cosalphasigma*(Xy-Sy))*(pp+1.0));

							if(value < minvaluelocal)
							{
								minvaluelocal = value;
								worstalphalocal = alpha;
								worstsigmalocal = sigma;
								worstpplocal = pp;
								worstchilocal = chi;
							}
						}
					}
				}
			}
		}
		
		#pragma omp critical
		{
			if(minvaluelocal < minvalue)
			{
				minvalue = minvaluelocal;
				worstalpha = worstalphalocal;
				worstsigma = worstsigmalocal;
				worstpp = worstpplocal;
				worstchi = worstchilocal;
			}
		}
	}
	
	cout << "Most problematic value: " << minvalue;
	cout << "\nAchieved at alpha = " << worstalpha << ", sigma = " << worstsigma << ", pp = " << worstpp << ", chi = " << worstchi;
	
	double endtime = omp_get_wtime();
	double cpu_time_used = ((double)(endtime - starttime));
	cout << "\n\nTime used: " << cpu_time_used << "s\n";
	
	ofstream outputfile;
	outputfile.open("output.txt");
	
	outputfile << "delta: " << delta << "\n\n";
	outputfile << "Most problematic value: " << minvalue;
	outputfile << "\nAchieved at alpha = " << worstalpha << ", sigma = " << worstsigma << ", pp = " << worstpp << ", chi = " << worstchi;
	outputfile << "\n\nTime used: " << cpu_time_used << "s\n";
	
	outputfile.close();

	return 0;
}