#include "dynamic_inverse_ising.hpp"

int DynamicInverseIsing::init()
{
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
    L = Ndim*Mdim;
    std::cout << "Ndim = " << Ndim << " Mdim = " << Mdim << "\n";

    output_file_name = "result.txt";
    input_path = "../evidence/";

    learned_s0 = new double[L];
    if(learned_s0==NULL)
    {
	std::cout<<"Allocating storage for learned_s0 FAILED!"<< "\n";
	return -1;
    }
    return 0;
}

int DynamicInverseIsing::getEvidence()
{
    std::string input_file_name;

    std::string in_file_list_name = input_path + "filelist.txt";
    FILE * in_file_list;
    in_file_list = fopen (in_file_list_name.c_str(), "r");
    if ( in_file_list == NULL ){
	std::cout << "open in_file_list error" << in_file_list_name << std::endl;
	getchar();
    }

    char in_data_file_name_part[200];
    std::string in_data_file_name;
    while(fscanf(in_file_list, "%s", in_data_file_name_part)!=EOF)
    {
	Evidence evidence;
	evidence.init(Ndim, Mdim);
	//std::cout << "here0.0" << "\n";

	in_data_file_name = input_path + in_data_file_name_part;
	std::cout << in_data_file_name << std::endl;

	const char * pszSrcFilename = in_data_file_name.c_str();
	GDALDataset * poSrcDS = (GDALDataset *) GDALOpen( pszSrcFilename, GA_ReadOnly );
	GDALRasterBand *poBand;

	poBand = poSrcDS->GetRasterBand(1);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  evidence.sigma, Ndim, Mdim, GDT_Int32, 0, 0 );
	poBand = poSrcDS->GetRasterBand(2);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  evidence.PM25Matrix, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poSrcDS->GetRasterBand(3);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  evidence.WindMatrix1, Ndim, Mdim, GDT_Float64, 0, 0 );
	poBand = poSrcDS->GetRasterBand(4);
	poBand->RasterIO( GF_Read, 0, 0, Ndim, Mdim,
			  evidence.WindMatrix2, Ndim, Mdim, GDT_Float64, 0, 0 );

	GDALClose( (GDALDatasetH) poSrcDS );
	for (int i=0; i<L; i++)
	{
	    if(evidence.PM25Matrix[i] < 0)
	    {
		evidence.valid_grid_list[i] = -1;
	    }
	    else
	    {
		evidence.valid_grid_list[i] = 1;
	    }
	}

	//std::cout << "here0" << "\n";
	evidence.get_pm25_mean();
	//std::cout << "here1" << "\n";
	//evidence.compute_ising_network();
	//std::cout << "here2" << "\n";
	evidence_list.push_back(evidence);
	//std::cout << "here3" << "\n";
	std::cout << "evidence list length: " << evidence_list.size() << "\n";
	//evidence.drop();
	//std::cout << "here4" << "\n";
    }
    std::cout << "evidence list length: " << evidence_list.size() << "\n";
    return 0;
}

DynamicInverseIsing::~DynamicInverseIsing()
{
}

int DynamicInverseIsing::update_using_x(double * x)
{
    for (int i=0; i<L; i++)
    {
	learned_s0[i] = x[i];
    }
    alpha = x[L];
    beta = x[L+1];
    lambda = x[L+2];
    return 0;
}

int DynamicInverseIsing::make_initial(double * x)
{
    return 0;
}

// int DynamicInverseIsing::pm25_to_sigma(double pm25)
// {
//     if (pm25 < 0) return -99999;
//     if (pm25 >= 200) return 1;
//     if (pm25 < 200) return -1;
//     std::cout << "error in fun" << "\n";
//     getchar();
// }

double DynamicInverseIsing::evaluate(const int N, const double * x, double * g)
{
    /* N is the number of unknowns to be solved. In our case N = L + 3 */
    double SL = 0.0;
    double A1, B1;
    double Pr;
    for (int i=0; i<N; i++)
    {
	g[i] = 0.0;
    }

    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	for(int r=0;r<L;r++) 
	{	
	    for ( int j=0; j<L; j++)
	    {
		iter->JJ[r][j] = -1.0;
	    }

	}
	for(int r=0;r<L;r++) 
	{	
	    iter->s[r] = x[r] + x[L+2]*iter->pm25_mean;
	    for ( int j=0; j<L; j++)
	    {
		if ( (r+1)%Mdim>0 )
		{
		    iter->JJ[r][j] = 0.01*x[L+1];
		}
		if ( r%Mdim>0 )
		{
		    iter->JJ[r][j] = 0.01*x[L+1];
		}
		if (r-Mdim >= 0)
		{
		    iter->JJ[r][j] = 0.01*x[L+1];
		}
		if (r+Mdim < L)
		{
		    iter->JJ[r][j] = 0.01*x[L+1];
		}
	    }
	}
	//std::cout << "hereEVA1" << "\n";
	for(int r=0;r<L;r++) 
	{	
	    if (iter->valid_grid_list[r] == 1)
	    {
		A1 = 0.0;
		B1 = 0.0;
		for ( int i=0; i<L; i++ )
		{
		    if ( (i==r+1) || (i==r-1) || (i==r-Mdim) || (i==r+Mdim))
		    {
			//printf("***>>>>>>*** \n");
			if ( iter->JJ[i][r] >= 0 )
			{
			    if (iter->JJ[r][i] < 0)
			    {
				printf("something wrong.\n");
				getchar();
			    }
			    A1 = A1 + iter->JJ[i][r] * iter->sigma[i] + 
				iter->JJ[r][i] * iter->sigma[i];
			    B1 = B1 + 0.01 * iter-> sigma[i] + 0.01 * iter->sigma[i];
			}
		    }
		}
		//printf("***---*** %.10f, %.10f\n", A1, iter->s[r]);
		A1 = -2.0*iter->sigma[r]*(iter->s[r] + A1);
		Pr = 1.0/(1.0+exp(A1));
		SL = SL - log(Pr);
		g[r] = g[r] - (1.0/Pr)*(exp(A1)*(2.0*iter->sigma[r]))/((1.0+exp(A1))*(1.0+exp(A1)));
		g[L] = 0.0;
		g[L+1] = g[L+1] - (1.0/Pr)*(exp(A1)*(2.0*iter->sigma[r]*B1))/((1.0+exp(A1))*(1.0+exp(A1)));
		g[L+2] = g[L+2] - (1.0/Pr)*(exp(A1)*(2.0*iter->sigma[r]*iter->pm25_mean))/((1.0+exp(A1))*(1.0+exp(A1)));
		//printf("***>>>*** %.10f, %.10f \n", A1, Pr);
	    }
	}
	//std::cout << "hereEVA10" << "\n";
    }
    SL = SL / evidence_list.size();
    for (int i=0; i<N; i++)
    {
	g[i] = g[i] / evidence_list.size();
    }    
    //printf("***>>>*** %.10f, %.10f \n", SL, g[0]);
    return SL;
}

int DynamicInverseIsing::output_result(const double * x)
{
    FILE * fp;
    if ( (fp = fopen(output_file_name.c_str(), "w")) == NULL )
    {
	std::cout << "file open failed. \n";
	getchar();
    }
    for (int i=0; i<L; i++)
    {
	learned_s0[i] = x[i];
	fprintf(fp, "%f\n", learned_s0[i]);
    }    
    alpha = x[L];
    fprintf(fp, "%f\n", alpha);
    beta = x[L+1];
    fprintf(fp, "%f\n", beta);
    lambda = x[L+2];
    fprintf(fp, "%f\n", lambda);
    fclose(fp);
    return 0;
}
