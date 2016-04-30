#include <iostream>
#include <dai/alldai.h>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"

using namespace dai;
using namespace std;

int pm25_to_sigma(double pm25)
{
    if (pm25 < 0) return -99999;
    if (pm25 >= 200) return 1;
    if (pm25 < 200) return -1;
    cout << "error in fun" << "\n";
    getchar();
}

FactorGraph img2fg( double * s, double ** JJ, int * sigma, int Ndim, int Mdim)
{
    vector<Var> vars;
    vector<Factor> factors;

    long int L = Ndim*Mdim;

    cout << "  Image width:  " << Ndim << endl;
    cout << "  Image height: " << Mdim << endl;
    // cout << "  Pairwise interaction strength:   " << JJ << endl;
    // Create a binary variable for each pixel
    for( long int i = 0; i < L; i++ )
        vars.push_back( Var( i, 2 ) );

    for( long int i = 0; i < L; i++ )
    {
	if ( sigma[i] == 1)
	{
	    if ( (i+1)%Mdim>0 )
	    {
		s[i+1] = s[i+1] + JJ[i][i+1]*sigma[i];
		JJ[i][i+1] = -99999.0;
		JJ[i+1][i] = -99999.0;
	    }
	    if ( i%Mdim>0 )
	    {
		s[i-1] = s[i-1] + JJ[i][i-1]*sigma[i];
		JJ[i][i-1] = -99999.0;
		JJ[i-1][i] = -99999.0;
	    }
	    if ( i-Mdim >= 0 )
	    {
		s[i-Mdim] = s[i-Mdim] + JJ[i][i-Mdim]*sigma[i];
		JJ[i][i-Mdim] = -99999.0;
		JJ[i-Mdim][i] = -99999.0;
	    }
	    if ( i+Mdim < L )
	    {
		s[i+Mdim] = s[i+Mdim] + JJ[i][i+Mdim]*sigma[i];
		JJ[i][i+Mdim] = -99999.0;
		JJ[i+Mdim][i] = -99999.0;
	    }
	}
	else if ( sigma[i] == -1)
	{
	    if ( (i+1)%Mdim>0 )
	    {
		s[i+1] = s[i+1] + JJ[i][i+1]*sigma[i];
		JJ[i][i+1] = -99999.0;
		JJ[i+1][i] = -99999.0;
	    }
	    if ( i%Mdim>0 )
	    {
		s[i-1] = s[i-1] + JJ[i][i-1]*sigma[i];
		JJ[i][i-1] = -99999.0;
		JJ[i-1][i] = -99999.0;
	    }
	    if (i-Mdim >= 0)
	    {
		s[i-Mdim] = s[i-Mdim] + JJ[i][i-Mdim]*sigma[i];
		JJ[i][i-Mdim] = -99999.0;
		JJ[i-Mdim][i] = -99999.0;
	    }
	    if (i+Mdim < L)
	    {
		s[i+Mdim] = s[i+Mdim] + JJ[i][i+Mdim]*sigma[i];
		JJ[i][i+Mdim] = -99999.0;
		JJ[i+Mdim][i] = -99999.0;
	    }
	}
	else
	{
	    if (sigma[i] != -99999)
	    {
		cout << "error sigma" << "\n";
		getchar();
	    }
	}
    }

    for( long int i = 0; i < L; i++ )
    {
	if ( sigma[i] == -99999)
	{
	    if ( (i+1)%Mdim>0 )
	    {
		if (JJ[i][i+1] > 0.0)
		{
		    factors.push_back( 
			createFactorIsing( vars[i], vars[i+1], JJ[i][i+1] ) );
		}
	    }
	    if ( i%Mdim>0 )
	    {
		if (JJ[i][i-1] > 0.0)
		{
		    factors.push_back( 
			createFactorIsing( vars[i], vars[i-1], JJ[i][i-1] ) );
		}
	    }
	    if (i-Mdim >= 0)
	    {
		if (JJ[i][i-Mdim] > 0.0)
		{
		    factors.push_back( 
			createFactorIsing( vars[i], vars[i-Mdim], JJ[i][i-Mdim] ) );
		}
	    }
	    if (i+Mdim < L)
	    {
		if (JJ[i][i+Mdim] > 0.0)
		{
		    factors.push_back( 
			createFactorIsing( vars[i], vars[i+Mdim], JJ[i][i+Mdim] ) );
		}
	    }
	    factors.push_back( createFactorIsing( vars[i], s[i] ) );
	}
    }
    // Create the factor graph out of the variables and factors
    cout << "Creating the factor graph..." << factors.size() << endl;
    return FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}

int main() 
{
    /* ***************** */
    /* Initialization    */
    /* ***************** */

    /* Get Ndim, Mdim */
    int Ndim, Mdim;
    GDALAllRegister();
    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    char **papszMetadata;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL )
	exit( 1 );
    papszMetadata = poDriver->GetMetadata();
    const char * pszSrcFilename = "SampleTiff.tif";
    GDALDataset *poSrcDS =
	(GDALDataset *) GDALOpen( pszSrcFilename, GA_ReadOnly );
    Ndim = poSrcDS->GetRasterXSize();
    Mdim = poSrcDS->GetRasterYSize();
    GDALClose( (GDALDatasetH) poSrcDS );
    long int L = Ndim*Mdim;
    cout << "Ndim = " << Ndim << " Mdim = " << Mdim << "\n";

    /* *********************** */
    /* Read network parameters */
    /* *********************** */
    double alpha, beta, lambda;
    double * s0;
    s0 = new double[L];
    if(s0==NULL)
    {
	std::cout<<"Allocating storage for s0 FAILED!"<< "\n";
	return -1;
    }
    std::string input_file_name = "result.txt";
    FILE * fp;
    if ( (fp = fopen(input_file_name.c_str(), "r")) == NULL )
    {
	std::cout << "file open failed. \n";
	getchar();
    }
    double temp;
    for ( int i=0; i<L; i++)
    {    
	fscanf(fp, "%lf\n", &temp);
	s0[i] = temp;
    }
    fscanf(fp, "%lf\n", &alpha);
    fscanf(fp, "%lf\n", &beta);
    fscanf(fp, "%lf\n", &lambda);
    fclose(fp);
    
    /* ***************** */
    /* Allocate Memory   */
    /* ***************** */
    double * PM25Matrix;
    double * WindMatrix1;
    double * WindMatrix2;
    PM25Matrix = new double[L];
    WindMatrix1 = new double[L];
    WindMatrix2 = new double[L];
    if(PM25Matrix==NULL | WindMatrix1==NULL | WindMatrix2==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    double ** JJ;
    double * s_total;
    int * sigma;
    s_total = new double[L];
    sigma = new int[L];
    if(s_total==NULL | sigma==NULL)
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

    /* ***************** */
    /* Data processing   */
    /* ***************** */
    std::string data_path = "../data/";
    std::string output_path = "../evidence/sigma_";
    std::string in_file_list_name = data_path + "filelist.txt";
    FILE * in_file_list;
    in_file_list = fopen (in_file_list_name.c_str(), "r");
    if ( in_file_list == NULL ){
      std::cout << "open in_file_list error" << in_file_list_name << std::endl;
	getchar();
    }

    char in_data_file_name_part[200];
    std::string in_data_file_name;
    /* For every file in the "filelist.txt", read in data and fill missing value */
    while(fscanf(in_file_list, "%s", in_data_file_name_part)!=EOF)
    {
	in_data_file_name = data_path + in_data_file_name_part;
	std::cout << in_data_file_name << std::endl;
	pszSrcFilename = in_data_file_name.c_str();
	poSrcDS = (GDALDataset *) GDALOpen( pszSrcFilename, GA_ReadOnly );
	GDALRasterBand *poBand;
	poBand = poSrcDS->GetRasterBand(2);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  PM25Matrix, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poSrcDS->GetRasterBand(3);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  WindMatrix1, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poSrcDS->GetRasterBand(4);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  WindMatrix2, Ndim, Mdim, GDT_Float64, 0, 0 );
    
	/* compute pm25_mean, which is a feature of the network */
	int nrecords = 0;
	double pm25_mean = 0.0;
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
    
	/* prepare network */
	for ( int i=0; i<L; i++)
	{
	    s_total[i] = 0.0 + lambda*pm25_mean;
	    sigma[i] = pm25_to_sigma(PM25Matrix[i]);
	    //sigma[i] = -99999;
	    //cout << i << " " << PM25Matrix[i] << " " << sigma[i] << "\n";
	    for ( int j=0; j<L; j++)
	    {
		JJ[i][j] = 0.01;
	    }
	}

	/* make factor graph, the outcome is a network of unknowns */
	FactorGraph fg = img2fg( s_total, JJ, sigma, Ndim, Mdim );
	//const char *file_fg = "FactorGraph.fg";
	// if( strlen( file_fg ) > 0 ) {
	//     cout << "Saving factor graph as " << file_fg << endl;
	//     fg.WriteToFile( file_fg );
	// }

	/* inference all the unknowns using Gibbs sampling */
	PropertySet gibbsProps;
	gibbsProps.set("maxiter", size_t(5000));   // number of Gibbs sampler iterations
	gibbsProps.set("burnin", size_t(4000));
	gibbsProps.set("verbose", size_t(0));
	Gibbs gibbsSampler( fg, gibbsProps );
	vector<size_t> mm;
	gibbsSampler.init();
	gibbsSampler.run();
	mm = gibbsSampler.state();
	//cout << "size mm" << mm.size() << "\n";

	/* generate the full state sigma: observed + inferred */
	int index = 0;
	for( size_t i = 0; i < L; i++ )
	{
	    if (sigma[i] == -99999)
	    {
		sigma[i] = 2*mm[i] - 1;
		//cout << i << " " << mm[i] << " " << sigma[i] << "\n";
	    }
	    //cout << "result " << i << " " << sigma[i] << "\n";
	}

	/* output full state */
	std::string out_file_name;
	out_file_name = output_path + in_data_file_name_part;	
	std::cout << out_file_name << std::endl;
	const char * pszDstFilename = out_file_name.c_str();
	GDALDataset *poDstDS = poDriver->CreateCopy( pszDstFilename, poSrcDS, FALSE,
					NULL, NULL, NULL );
	poBand = poDstDS->GetRasterBand(1);
	//poBand->RasterIO( GF_Write, 0, 0, Ndim, Mdim,
	//PM25Matrix, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand->RasterIO( GF_Write, 0, 0, Ndim, Mdim,
			  sigma, Ndim, Mdim, GDT_Int32, 0, 0 );
	poBand = poDstDS->GetRasterBand(2);
	poBand->RasterIO( GF_Write, 0, 0, Ndim, Mdim,
			  PM25Matrix, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poDstDS->GetRasterBand(3);
	poBand->RasterIO( GF_Write, 0, 0, Ndim, Mdim,
			  WindMatrix1, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poDstDS->GetRasterBand(4);
	poBand->RasterIO( GF_Write, 0, 0, Ndim, Mdim,
			  WindMatrix2, Ndim, Mdim, GDT_Float64, 0, 0 );
	GDALClose( (GDALDatasetH) poDstDS );
	GDALClose( (GDALDatasetH) poSrcDS );
    }
    return 0;
}
