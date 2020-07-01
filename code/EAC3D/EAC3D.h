#ifndef __EAC3D_DOT_H__
#define __EAC3D_DOT_H__




#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <vector>





#include "SuperLU\SRC\slu_ddefs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>


#define TID 0   //Temperature id
#define HID 1   //Water content id
#define rD 2   //degree of reaction id
#define STRAINID 3   //strain id
#define STRESSID 4   //stress id


#define BND_DIR 6   //boundary condition DIRICHLET 
#define BND_NEU 7   //boundary condition NEUMANN 
#define BND_CON 8   //boundary condition CONNECTION 
enum ThermalDiffModel {scalar};

typedef double(*func_ptr_dd)(double);

class  EAC3D{
    
private:
    double x0, y0, z0;
    double Lx, Ly, Lz;
    int nx, ny, nz;
    double dx, dy, dz;
	std::vector<double>*** Strain;
	std::vector<double>*** Stress;
	std::vector<double> TimeVector;
    double**** F;
	double**** Finit;
	double**** Fstar;
    double**** RHSnm1;
    
    double**** RHS;

	bool isFirstStep;
	int innerRK2Step;
    double time; 
    
    double dt;
	int DIM_ARRAY;

	func_ptr_dd *DCoef;

    //mass density [kg/m^3]
    double rho;
    //specific heat capacity [J/(kg.K)]
    double cp;
    
    
    double (*bnd_T_xs)(double, double,double, double);
    double (*bnd_T_xf)(double, double,double, double);
    double (*bnd_T_ys)(double, double,double, double);
    double (*bnd_T_yf)(double, double,double, double);
    double (*bnd_T_zs)(double, double,double, double);
    double (*bnd_T_zf)(double, double,double, double);
	func_ptr_dd D_T;

    int type_bnd_T_xs;
    int type_bnd_T_xf;
    int type_bnd_T_ys;
    int type_bnd_T_yf;
    int type_bnd_T_zs;
    int type_bnd_T_zf;
    
    
    double (*bnd_H_xs)(double, double,double, double);
    double (*bnd_H_xf)(double, double,double, double);
    double (*bnd_H_ys)(double, double,double, double);
    double (*bnd_H_yf)(double, double,double, double);
    double (*bnd_H_zs)(double, double,double, double);
    double (*bnd_H_zf)(double, double,double, double);

	double (*DStrainODH)(double);
	double (*DStrainODT)(double);
	double (*J_relax)(double, double);
	double (*Restrain_func)(double, double, double);

	func_ptr_dd D_H;

    int type_bnd_H_xs;
    int type_bnd_H_xf;
    int type_bnd_H_ys;
    int type_bnd_H_yf;
    int type_bnd_H_zs;
    int type_bnd_H_zf;

    
    void compute_RHS();
    
    
    void set_BNDXS(double (*bnd)(double, double , double, double), int id, int type);
    void set_BNDXF(double (*bnd)(double, double , double, double), int id, int type);
    
    void set_BNDYS(double (*bnd)(double, double , double, double), int id, int type);
    void set_BNDYF(double (*bnd)(double, double , double, double), int id, int type);
    
    void set_BNDZS(double (*bnd)(double, double , double, double), int id, int type);
    void set_BNDZF(double (*bnd)(double, double , double, double), int id, int type);
    
    
public :
    EAC3D(double Lx, double Ly, double Lz, int nx, int ny, int nz);
    ~EAC3D();
    
    void  set_T_BNDXS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_T_BNDXF(double (*bnd)(double, double , double, double), int bnd_type);
    
    void  set_T_BNDYS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_T_BNDYF(double (*bnd)(double, double , double, double), int bnd_type);
    
    void  set_T_BNDZS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_T_BNDZF(double (*bnd)(double, double , double, double), int bnd_type);
    
    void  set_H_BNDXS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_H_BNDXF(double (*bnd)(double, double , double, double), int bnd_type);
    
    void  set_H_BNDYS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_H_BNDYF(double (*bnd)(double, double , double, double), int bnd_type);
    
    void  set_H_BNDZS(double (*bnd)(double, double , double, double), int bnd_type);
    void  set_H_BNDZF(double (*bnd)(double, double , double, double), int bnd_type);

	void set_DStrainODH(double (*func)(double));
	void set_DStrainODT(double (*func)(double));
    
	void set_Restrain(double(*func_restrain)(double, double, double));
	void set_Relaxation(double(*func_J_relax)(double, double));
    void  set_T_lambda(double(*DiffusionCoef)(double));
	void  set_H_lambda(double(*DiffusionCoef)(double));
    void  set_rho(double val);
    void  set_cp(double val);
    
    void  saveData(std::string s, int bnd_id);
    void  saveMesh();
    
    double set_dt(double dt);
    void  timeMarching(bool onlyRK2);
    
    void  init_field(double (*init)(double, double,double), int bnd_id);
	void set_NBEquations(int NBEQ);
    
	void implicitTimeMarching(double CN_factor, int nbIter);
	void computeStrainStress();
    
    
};

#endif
