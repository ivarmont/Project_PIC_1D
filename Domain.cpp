#include "Domain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Functions.h"

Domain::Domain(){
	Nx=1;
	dx=0.0;
	x0=0.0;
	xl= 0.0;
	x_s0=0.0;	
	x_s1=0.0;

	rho = new double[Nx];
	phi = new double[Nx];
	ef = new double[Nx];

	res = new double[Nx];
	cg_dir = new double[Nx];
	Acg_dir = new double[Nx];

	};
Domain::Domain( int N_grid_points, double cell_size){
	Nx = N_grid_points;
	dx= cell_size;

	x0=0.0;
	xl = (N_grid_points-1)*dx;
	x_s0 = x0 + 0.01*xl;
	x_s1 = x0 + 0.25*xl;		
	/*
	n_e = new double[Nx];
	n_pi = new double[Nx];
	n_ni = new double[Nx];
	*/
	rho = new double[Nx];
	phi = new double[Nx];
	ef = new double[Nx];

	res = new double[Nx];
	cg_dir = new double[Nx];
	Acg_dir = new double[Nx];
	for (int i=0;i<Nx;i++)
	{
		res[i]= 0.0;
		cg_dir[i]= 0.0;
		Acg_dir[i]= 0.0;	
	}

	};
	
Domain::~Domain(){
	
	};

double Domain::get_domain_length(void){
	return xl;

};

double Domain::get_domain_start(void){
	return x0;

};

double Domain::get_domain_end(void){
	return (x0 + xl);

};

double Domain::get_domain_source_begin(void){
	return x_s0;

};

double Domain::get_domain_source_end(void){
	return x_s1;

};


void Domain::clean_charge(void){
	 
	for(int i =0;i<Nx;i++){
		rho[i] = 0.0;
	}
};

void Domain::set_charge(double position, double charge){
	double ld=(position-x0)/dx;
	int li=ld;
	double distance = ld-li; 

	rho[li]+=charge*(1.0-distance)/dx;
	rho[li + 1]+=charge*(distance)/dx;
};

void Domain::set_charge_2(double position, double charge){
	double ld=(position-x0)/dx;
	int li=ld;
	double distance = ld-li; 
	double vol=1.0/(4.0*dx);
	
	if(li>0 && li<(Nx-2)){ 
		rho[li-1] += vol*charge*(1.0-distance);
		rho[li] += vol*charge*(2.0-distance);
		rho[li+1] += vol*charge*(1.0 + distance);
		rho[li+2] += vol*charge*(distance);
	}else{
		rho[li] += vol*charge*(1-distance);
		rho[li+1] += vol*charge*(distance);	
	}
};



bool Domain::solve_potential_GS(double phi_0, double phi_N){
	double L2;
	double dx2=dx*dx;
	double EPS0= 8.85418782e-12;

	/*Initialize boundaries*/
	phi[0]=phi_0;
	phi[Nx-1]=phi_N;
	
	for(int solver_it=0;solver_it<1000000;solver_it++)
	{
		for(int i=1;i<Nx-1;i++)
		{
			/*SOR*/	
			double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]/EPS0 );
			phi[i] = phi[i] + 1.4*(g-phi[i]);
		}
	
	/*check for convergence*/
	if(solver_it%25==0)
	{
		double sum=0;
		for(int i=1;i<Nx-1;i++)
		{
			double R = -rho[i]/EPS0 - (phi[i-1] - 2.0*phi[i] + phi[i+1])/dx2;
			sum+=R*R;
		}
		L2 =sqrt(sum)/Nx;
		if(L2<1e-10)
		{
			return true;
		}	

	}

	}

	std::cout<<"Gauss-Seidel solver failed to converge, L2 = " << L2 <<std::endl;
	return false;
};

bool Domain::solve_potential_CG( double phi_0, double phi_N)
{
	double error;
	double dx2 = dx*dx;	/*precompute*/
	double tolerance=1e-18;
	double r2=0.0;
	double r2_new=0.0;	
	double pAp=0.0;
	double alpha=0.0;
	double beta=0.0;
	double EPS_0= 8.85418782e-12;
	
	/*initialize boundaries*/
	phi[0]=phi_0;	
	phi[Nx-1]=phi_N;
	
	for(int i=1;i<(Nx-1);i++)
	{
		res[i]= -(rho[i]/EPS_0)  - (phi[i-1] - 2.0*phi[i] + phi[i+1] )/dx2;
		cg_dir[i]=res[i];	
	}
	
	if(dot_product(res,res, Nx)>0){
	error=0.0;

	/*solve potential*/
	for (int solver_it=0;solver_it<4000;solver_it++)
	{
		for(int i=1;i<(Nx-1);i++)
		{
			Acg_dir[i]=(cg_dir[i-1] - 2.0*cg_dir[i] + cg_dir[i+1] )/dx2;				
		}
	
	
		r2= dot_product(res,res, Nx);
		pAp=dot_product(cg_dir,Acg_dir, Nx);
		alpha=r2/pAp;
		
		for(int i=1;i<(Nx-1);i++)
		{
			phi[i]=phi[i] + alpha*cg_dir[i];
			res[i]=res[i]-alpha*Acg_dir[i];	
		}
	
		r2_new= dot_product(res,res, Nx);
		beta=r2_new/r2;
		for(int i=1;i<(Nx-1);i++)
		{
			cg_dir[i]=res[i] + beta*cg_dir[i];	
		}	


			error=sqrt(r2_new)/Nx;
			if (error<tolerance) {		
				return true;
			}
	
	}
	std::cout<< "CG solver failed to converge, error= "<< error <<std::endl;
	return false;	
	}
	return true;

}







void Domain::compute_ef(){
	for( int i =1; i <Nx-1; i++) ef[i]= -(phi[i+1]-phi[i-1])/(2.0*dx);
	ef[0] = -(phi[1] - phi[0])/dx;
	ef[Nx-1]= -( phi[Nx-1] - phi[Nx-2])/dx; 
};


double Domain::interpolate_ef(double xpos){
	double ld=(xpos-x0)/dx;
	int li=ld;
	double distance = ld-li;
	double ef_inter=0.0;
	
	ef_inter= ef[li]*(1.0-distance) +  ef[li + 1]*(distance);
	return ef_inter;	
	
	 
};

double Domain::interpolate_ef_2(double xpos){
	double ld=(xpos-x0)/dx;
	int li=ld;
	double distance = ld-li;
	double ef_inter=0.0;
	
	
	if(li>0 && li<(Nx-2)){ 
		ef_inter = (ef[li-1]*(1.0-distance) + ef[li]*(2.0-distance) + ef[li+1]*(1.0+distance) + ef[li+2]*(distance))*0.25; 
	}else{
		ef_inter= ef[li]*(1.0-distance) + ef[li+1]*(distance);
	}

	return ef_inter;	

	 
};

double * Domain::get_potential(){

	return phi;

}


void Domain::print_charge_density(){
	std::string filename="charge.dat";
	std::fstream file;
	file.open(filename, std::ios::out);	
	for (int i=0;i<Nx;i++){	
		file<< i*dx <<"\t"<< rho[i] << std::endl; 
	}
	file.close();

};



void Domain::print_potential(){
	std::string filename="potential.dat";
	std::fstream file;
	file.open(filename, std::ios::out);	
	for (int i=0;i<Nx;i++){	
		file<< i*dx <<"\t"<< phi[i] << std::endl; 
	}
	file.close();

};


void Domain::print_ef(){
	std::string filename="ef.dat";
	std::fstream file;
	file.open(filename, std::ios::out);	
	for (int i=0;i<Nx;i++){	
		file<< i*dx <<"\t"<< ef[i] << std::endl; 
	}
	file.close();

};

	
