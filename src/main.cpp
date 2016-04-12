//=============================================================================
//  iSpectra
//=============================================================================
//
//  Copyright 2012 Chun Shen & Zhi Qiu
//    shen.201@asc.ohio-state.edu
//    qiu.24@asc.ohio-state.edu
//
//        Date: 03/2012

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

#include "./main.h"
#include "./Table.h"
#include "./readindata.h"
#include "./emissionfunction.h"
#include "./arsenal.h"
#include "./ParameterReader.h"

using namespace std;

int main(int argc, char *argv[]) {
    cout << endl
         << "                      iSpectra                   " << endl
         << endl
         << "  Ver 1.2   ----- Zhi Qiu & Chun Shen, 03/2012   " << endl;
    cout << endl
         << "**********************************************************"
         << endl;
    display_logo(2);  // Hail to the king~
    cout << endl << "**********************************************************"
         << endl << endl;

    // Read-in parameters
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->readFromArguments(argc, argv);
    paraRdr->echo();

    // Chun's input reading process
    string path = "results";

    // load freeze out information
    read_FOdata freeze_out_data(paraRdr, path);

    int FO_length = 0;
    FO_length = freeze_out_data.get_number_of_freezeout_cells();
    cout << "total number of cells: " <<  FO_length << endl;

    FO_surf* FOsurf_ptr = new FO_surf[FO_length];
    for (int i = 0; i < FO_length; i++)
        for (int j = 0; j < Maxparticle; j++)
            FOsurf_ptr[i].particle_mu[j] = 0.0e0;

    freeze_out_data.read_in_freeze_out_data(FO_length, FOsurf_ptr);

    // read the chemical potential on the freeze out surface
    particle_info *particle = new particle_info[Maxparticle];
    int Nparticle = freeze_out_data.read_in_chemical_potentials(
                                        path, FO_length, FOsurf_ptr, particle);

    cout << endl << " -- Read in data finished!" << endl << endl;

    // check chemical potentials
    // for(int j=0; j<Nparticle; j++)
    // cout << FOsurf_ptr[0].Edec << "   " << FOsurf_ptr[0].particle_mu[j]
    //      << endl;
    // exit(0);

    // Next, Zhi's turn...
    // skip others except for these particle
    Table chosen_particles("EOS/chosen_particles.dat");
    Table pT_tab("tables/pT_gauss_table.dat");  // pt position and weight table
    // phi position and weight table
    Table phi_tab("tables/phi_gauss_table.dat");
    // eta uniform dist table
    Table eta_tab("tables/eta_gauss_table_15.dat");
    EmissionFunctionArray efa(paraRdr, 0.0, &chosen_particles, &pT_tab,
                              &phi_tab, &eta_tab, particle, Nparticle,
                              FOsurf_ptr, FO_length);

    efa.calculate_dN_ptdptdphidy_and_flows_4all(9);

    delete [] FOsurf_ptr;
    delete paraRdr;
}
