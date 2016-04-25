#ifndef EVIDENCE_HPP_
#define EVIDENCE HPP_

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

class Evidence
{
public:
    double ** JJ;
    double * s;
    int * sigma;
    double * PM25Matrix;
    double * WindMatrix1;
    double * WindMatrix2;
    double pm25_mean;
    int Ndim, Mdim;
    long int L;

public:
    int drop();
    int init(int N, int M);
    int compute_ising_network();
    int get_pm25_mean();
};
#endif /* EVIDENCE_HPP_ */
