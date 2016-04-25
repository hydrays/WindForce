#include "evidence.hpp"

int Evidence::drop()
{
    delete s;
    for (int i=0; i<L; i++)
    {
	delete JJ[i];
    }
    delete JJ;
    delete sigma;
    delete PM25Matrix;
    delete WindMatrix1;
    delete WindMatrix2;
};

int Evidence::init(int N, int M)
{
    Ndim = N;
    Mdim = M;
    L = Ndim * Mdim;
	    
    sigma = new int[L];
    if(sigma==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    PM25Matrix = new double[L];
    WindMatrix1 = new double[L];
    WindMatrix2 = new double[L];
    if(PM25Matrix==NULL | WindMatrix1==NULL | WindMatrix2==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    s = new double[L];
    if(s==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    JJ = new double*[L];
    if(JJ==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    for ( int i=0; i<L; i++)
    {
	JJ[i] = new double[L];
	if(JJ[i]==NULL)
	{
	    std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	    return -1;
	}
    }
};

int Evidence::compute_ising_network()
{
    for ( int i=0; i<L; i++)
    {
	s[i] = 0.0;
	//ss[i] = s[i];
	//sigma[i] = pm25_to_sigma(PM25Matrix[i]);
	//sigma[i] = -99999;
	//std::cout << i << " " << PM25Matrix[i] << " " << sigma[i] << "\n";
	for ( int j=0; j<L; j++)
	{
	    JJ[i][j] = 0.01;
	}
    }
    return 0;
};

int Evidence::get_pm25_mean()
{
    int nrecords;
    pm25_mean = 0.0;
    nrecords = 0;
    for ( int i=0; i<L; i++)
    {
	if (PM25Matrix[i] > 0.0)
	{
	    pm25_mean = pm25_mean + PM25Matrix[i];
	    nrecords = nrecords + 1;
	}
    }
    pm25_mean = pm25_mean / nrecords;
    pm25_mean = 0.01*pm25_mean;
    return 0;
};
