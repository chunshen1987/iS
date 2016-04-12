// Copyright 2012 Chun Shen and Zhi Qiu
#ifndef SRC_EMISSIONFUNCTION_H_
#define SRC_EMISSIONFUNCTION_H_

#include <string>
#include <vector>
#include "./Table.h"
#include "./main.h"
#include "./ParameterReader.h"

using namespace std;

class EmissionFunctionArray {
 private:
    ParameterReader* paraRdr;

    int CALCULATEDED3P;
    int INCLUDE_BULKDELTAF, INCLUDE_DELTAF;
    int INCLUDE_MUB;
    int bulk_deltaf_kind;
    int GROUPING_PARTICLES;
    double PARTICLE_DIFF_TOLERANCE;
    int USE_HISTORIC_FORMAT;
    int F0_IS_NOT_SMALL;

    int hydro_mode;
    double particle_y;
    Table *pT_tab, *phi_tab, *eta_tab;
    int pT_tab_length, phi_tab_length, eta_tab_length;
    long FO_length;
    Table *dN_ptdptdphidy;
    Table *dE_ptdptdphidy;
    int number_of_chosen_particles;

    // has length Nparticle, 0 means miss, 1 means include
    int *chosen_particles_01_table;

    // store particle index;
    // the sampling process follows the order specified by this table
    int *chosen_particles_sampling_table;

    int Nparticles;
    particle_info* particles;
    FO_surf* FOsurf_ptr;
    // store the last particle index being used by calculate_dN_ptdptdphidy()
    int last_particle_idx;
    bool particles_are_the_same(int, int);

    // array for bulk delta f coefficients
    Table *bulkdf_coeff;

 public:
    EmissionFunctionArray(ParameterReader* paraRdr_in, double particle_y_in,
                          Table* chosen_particle, Table* pT_tab_in,
                          Table* phi_tab_in, Table* eta_tab_in,
                          particle_info* particles_in, int Nparticles,
                          FO_surf* FOsurf_ptr_in, long FO_length_in);
    ~EmissionFunctionArray();

    void calculate_dN_ptdptdphidy(int particle_idx);
    void calculate_dN_ptdptdphidy_3D(int particle_idx);
    void write_dN_ptdptdphidy_toFile();
    string dN_ptdptdphidy_filename;  // where to save
    string dE_ptdptdphidy_filename;  // where to save

    void calculate_flows(int to_order, string, string);
    void calculate_Energyflows(int to_order, string, string);

    string flow_differential_filename_old, flow_integrated_filename_old;
    string flow_differential_filename, flow_integrated_filename;
    string energyflow_differential_filename_old;
    string energyflow_integrated_filename_old;
    string energyflow_differential_filename, energyflow_integrated_filename;

    void calculate_dN_ptdptdphidy_and_flows_4all(int to_order = 9);
    void getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients);
};

#endif  // SRC_EMISSIONFUNCTION_H_
