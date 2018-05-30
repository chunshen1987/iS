// Copyright 2012 Chun Shen and Zhi Qiu
#ifndef SRC_READINDATA_H_
#define SRC_READINDATA_H_

#include <fstream>
#include <string>

#include "./main.h"
#include "./ParameterReader.h"

using namespace std;

typedef struct {
    int monval;     // Monte Carlo number according PDG
    string name;
    double mass;
    double width;
    int gspin;      // spin degeneracy
    int baryon;
    int strange;
    int charm;
    int bottom;
    int gisospin;   // isospin degeneracy
    int charge;
    int decays;     // amount of decays listed for this resonance
    int stable;     // defines whether this particle is considered as stable
    int decays_Npart[Maxdecaychannel];
    double decays_branchratio[Maxdecaychannel];
    int decays_part[Maxdecaychannel][Maxdecaypart];
    int sign;       // Bose-Einstein or Dirac-Fermi statistics
} particle_info;

typedef struct {
    double tau, xpt, ypt, eta;
    double cosh_eta, sinh_eta;
    double da0, da1, da2, da3;
    double u0, u1, u2, u3;
    double Edec, Tdec, Pdec;
    double Bn, muB, muS;
    double pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
    double bulkPi;
    double particle_mu[Maxparticle];
} FO_surf;

class read_FOdata {
 private:
    ParameterReader* paraRdr;
    string path;
    int mode;
    int turn_on_bulk;
    int turn_on_muB;
    int n_eta_skip;
    int IEOS_music;
    bool surface_in_binary;

 public:
    read_FOdata(ParameterReader* paraRdr_in, string path);
    ~read_FOdata();

    int get_number_of_freezeout_cells();
    int get_number_of_lines_of_binary_surface_file(string filename);
    void read_in_freeze_out_data(int length, FO_surf* surf_ptr);
    int read_in_chemical_potentials(
        string path, int FO_length, FO_surf* surf_ptr,
        particle_info* particle_ptr);
    void read_decdat(int length, FO_surf* surf_ptr);
    void read_surfdat(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_VISH2p1(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_MUSIC(int length, FO_surf* surf_ptr);
    void read_FOsurfdat_MUSIC_boost_invariant(int length,
                                              FO_surf* surf_ptr);
    void read_FOsurfdat_hydro_analysis_boost_invariant(int length,
                                                       FO_surf* surf_ptr);
    void read_decdat_mu(int FO_length, int N_stable, double** particle_mu);
    void read_chemical_potentials_music(int FO_length, FO_surf* FOsurf_ptr,
                                        int N_stable, double** particle_mu);
    int read_resonances_list(particle_info* particle);
    void calculate_particle_mu(int Nparticle, FO_surf* FOsurf_ptr,
                               int FO_length, particle_info* particle,
                               double** particle_mu);
    void regulate_surface_cells(int length, FO_surf* surf_ptr);
    void regulate_Wmunu(double u[4], double Wmunu[4][4],
                        double Wmunu_regulated[4][4]);
};

#endif  // SRC_READINDATA_H_
