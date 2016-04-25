#ifndef DYNAMIC_INVERSE_ISING_HPP_
#define DYNAMIC_INVERSE_iSING HPP_

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "evidence.hpp"

class DynamicInverseIsing
{
public:
    double lambda;
    double beta;
    int Ndim, Mdim;
    long int L;
    double * learned_s;
    std::string output_file_name;
    std::string input_path;
    std::vector<Evidence> evidence_list;

public:
    int init();
    ~DynamicInverseIsing();
    int getEvidence();
    int update_using_x(double *);
    //int compute_ising_network();
    int make_initial(double *);
    //int pm25_to_sigma(double);
    //int get_pm25_mean();
    double evaluate(const int, const double *, double *);
    int output_result();
};

#endif /* DYANMIC_INVERSE_ISING_HPP_ */
