#include <iostream>
#include <fstream>	
#include "Particle.h"
#include "Domain.h"
#include "Functions.h"
#include <cmath>
#include <vector>

struct Specie
{
double charge;
double mass;
double weight;
double temperature;
double density;		


int Np; //Number of particles
int Np_alloc;
int Np_lost;

Particle *part;

double *n_avg;

};



//Add particle function
void Add_Particle(Specie *specie, double x, double v)
{

	if( specie->Np > specie->Np_alloc )
	{
		std::cout<< "Too many particles" << std::endl; 
	}
	else{
		//I add the particle in the last position
		specie->part[specie->Np].set_position(x);
		specie->part[specie->Np].set_velocity(v);
		specie->part[specie->Np].set_is_alive(true);	
		//I increment the counter of particles
		specie->Np++;
	}
}

void ScatterSpecie( Specie *specie, Domain *dom){
	for(int i=0;i < specie->Np; i++)	
		if (specie->part[i].is_alive() ) dom->set_charge( specie->part[i].get_position() ,  specie->weight*specie->charge);	
} 

void ScatterSpecie_2( Specie *specie, Domain *dom){
	for(int i=0;i < specie->Np; i++)	
		if (specie->part[i].is_alive() ) dom->set_charge_2( specie->part[i].get_position() ,  specie->weight*specie->charge);	
} 

void RewindSpecies( Specie *specie, Domain *dom, double dt){

	double ef=0.0;
	double qm=0.0;
	double force=0.0;	

	qm=specie->charge/specie->mass;
	specie->Np_lost=0;

	for(int i=0;i<specie->Np;i++){			
		if (specie->part[i].is_alive() ){		
			ef = dom->interpolate_ef_2( specie->part[i].get_position() );	
			force =-0.5*ef*qm;			
			specie->part[i].move_particle_LF(force, dt);			
		}

	}

}


void Push_specie( Specie *specie, Domain *dom, double dt){

	double ef=0.0;
	double qm=0.0;
	double force=0.0;	

	qm=specie->charge/specie->mass;
	specie->Np_lost=0;

	for(int i=0;i<specie->Np;i++){			
		if (specie->part[i].is_alive() ){		
			ef = dom->interpolate_ef_2( specie->part[i].get_position() );	
			force =ef*qm;			
			specie->part[i].move_particle_LF(force, dt);		
			
			if ( (specie->part[i].get_position() <  dom->get_domain_start()) || (specie->part[i].get_position() >  dom->get_domain_end())    ){
				specie->part[i].set_is_alive( false);		
				specie->Np_lost++;		
			}	
		}

	}



}


void thermalization(Specie *specie, Domain *dom ,double v_th, double dt)
{
	double v=0.0;
	double aux_rand=0.0;
	double w_e_reset=1.0e9;
	double prob_of_coll=1.0-exp(-w_e_reset*dt);

	
	for(int i=0;i<specie->Np;i++)
	{

		if( (specie->part[i].get_position() > dom->get_domain_source_begin() ) && (specie->part[i].get_position() < dom->get_domain_source_end() )   ){
			aux_rand=rnd();				
				if( aux_rand < prob_of_coll )
				{
					v = SampleVel2(v_th);
					specie->part[i].set_velocity(v);
				}		
		}
		
	}
}	


void Reorder_Specie(Specie *specie)
{	
	Particle *part_alive=new Particle[specie->Np];
	//Loop to move over each particlea
	int p=0;
	int m=0;	
	for(p=0 ; p < specie->Np; p++)
	{
		if(specie->part[p].is_alive())
		{
			part_alive[m]=specie->part[p];
			m++;
		}
	}

	specie->Np=m;
	p=0;
	for( p=0 ; p < specie->Np; p++)
	{
		specie->part[p]=part_alive[p];	
	}
	delete [] part_alive;
}


void reinjection_B(Specie *electrons, Specie *ions, Domain *dom ,double v_the, double v_thi)
{
	double x=0.0;
	double v=0.0;
	double xs0=dom->get_domain_source_begin();
	double xs1=dom->get_domain_source_end();
	for(int i=0;i<ions->Np_lost;i++)
	{
			x = xs0 + rnd()*(xs1-xs0);//p*delta_ions;//
			v = SampleVel2(v_thi);
			Add_Particle(ions, x,v);
	
			v = SampleVel2(v_the);
			Add_Particle(electrons, x,v);
	}
}




void save_Specie(Specie *specie, int i)
{	
	std::string filename="Specie_"+ std::to_string(i) +".dat";
	std::fstream file;
	file.open(filename, std::ios::out);	
	

	file<< specie->charge << std::endl; 
	file<< specie->mass << std::endl; 
	file<< specie->weight << std::endl;
	file<< specie->temperature << std::endl;
	file<< specie->density << std::endl;	
	file<< specie->Np << std::endl;
	file<< specie->Np_alloc << std::endl;

	file.close();


	filename="particles_spc_" + std::to_string(i) +".dat";
	file.open(filename, std::ios::out);

	
	
	//Loop to move over each particle
	int p=0;
		
	for(p=0 ; p < specie->Np; p++)
	{	
		file<<	specie->part[p].get_position() << "\t" << specie->part[p].get_velocity()<< "\t" << specie->part[p].is_alive()<<std::endl;
	}			
	file.close();
}

void load_Specie(Specie *specie, int i)
{	
	double x,v;
	bool state_alive;
	std::vector<double> line_d(5,0);
	std::vector<int> line_i(2,0);	
	std::string filename="Specie_"+ std::to_string(i) +".dat";
	
	
	std::fstream file;
	file.open(filename, std::ios::in);	
	
	if( file.is_open()){
		std::cout << filename<< " opened " << std::endl;

		while(!file.eof()){
			for (int i=0;i<7;i++)
			{
				if(i<5)
			 		file>>line_d[i];		
				else 
					file>>line_i[ i-5 ];
			}
			specie->charge = line_d[0]; 
			specie->mass = line_d[1]; 
			specie->weight = line_d[2];
			specie->temperature = line_d[3];			
			specie->density = line_d[4];
			specie->Np = line_i[0];
			specie->Np_alloc = line_i[1];
			std::cout<< " sp_charge = " <<specie->charge << " Np = " << specie->Np << " Np_alloc  " <<  specie->Np_alloc <<std::endl; 	
		}
		std::cout << " Specie data read" << std::endl;
	//	file.close();
		}
/*
	filename="particles_spc_" + std::to_string(i) +".dat";
	file.open(filename, std::ios::in);

	if( file.is_open()){
		std::cout << filename<< " opened " << std::endl;

		specie->part = new Particle[2*specie->Np_alloc];
		//Loop to move over each particle
		int p=0;
		x=0.0;
		v=0.0;
	
		while(!file.eof()){
			//if(p <= specie->Np){
				file >>	x >> v >> state_alive;
				specie->part[p].set_position(x); 
				specie->part[p].set_velocity(v)	;
				specie->part[p].set_is_alive(true);
				p++;
			//}
		}	
		
		//if(p!=specie->Np) std::cout << "Error loading particles -- " << "Np = "<< specie->Np<< "  Particles loaded = "<< p << std::endl;	
			
		file.close();
	}
*/
}




double XtoL(double pos, double x0, double dx)
{	
	double li=(pos-x0)/dx;
	return li;
}


void scatter2(double lc, double value, double *field, int Nx, double vol){
	int i = (int)lc;
	double di= lc-i;
	
	if(i>0 && i<(Nx-2)){ 
		field[i-1] += vol*value*(1.0-di);
		field[i] += vol*value*(2.0-di);
		field[i+1] += vol*value*(1.0 + di);
		field[i+2] += vol*value*(di);
	}else{
		field[i] += vol*value*(1-di);
		field[i+1] += vol*value*(di);	
	}
}


void ScatterSpecies_den(Specie *species, double Nx ,double dx, double x0){

	double vol=1.0/(4.0*dx);

	//Proyecto las partículas sobre la grilla
	for(int p=0;p<species->Np;p++){
		double lc=XtoL(species->part[p].get_position(), x0, dx);
		scatter2(lc , species->weight, species->n_avg, Nx, vol);
	}

}



int main(void){	
	
	static const double Kb=1.38065e-23; // J/K, Boltzmann constant
	static const double EV_TO_K=11604.52; // 1eV in Kelvin, QE/K	
	static const double QE =1.602176565e-19;	 // C, elementary charge
	static const double AMU=1.660538921e-27; // kg, atomic mass unit	
	static const double ME =9.10938215e-31;		// kg, electron mass 
	
	static const int NUM_TS=100000;

	//Domain parameters	
	int Nx=101;
	double dx=5.0e-5;
	double dt=2.0e-11;
	double density=1.0e16;
	double weight=5e9;


	bool re_start=false;
	double x=0.0;
	double v=0.0;	


	int Np_initial=0;
	double *phi_avg = new double[Nx];
	double *phi_aux = new double[Nx];
	std::string conver_f="convergence.dat";
	std::fstream file_c;
	file_c.open(conver_f, std::ios::out);	
	
	std::string average_f="average.dat";
	std::fstream file_avg;
	file_avg.open(average_f, std::ios::out);

	//Create domain
	Domain sim = Domain(Nx, dx);
	std::cout<< "Domain length = " << sim.get_domain_length() << std::endl;

	//Create species

	Specie electrons;
	Specie ions;	

	re_start=true;
	if(re_start == false)
	{
		electrons.charge= -QE; 
		electrons.mass=ME ;
		electrons.weight= weight ;
		electrons.temperature= 1.0;
		electrons.density= density ;
		Np_initial=(Nx-1)*dx*electrons.density/electrons.weight;	
		std::cout<<" Initial amount of e Particles simulated " << Np_initial <<std::endl;
		electrons.Np=0;
		electrons.Np_alloc=2*Np_initial;
		electrons.part = new Particle[2*electrons.Np_alloc];
		electrons.n_avg = new double[Nx];		

	
		ions.charge= QE ; 
		ions.mass= 1.0*AMU ;
		ions.weight= weight ;
		ions.temperature = 1.0;
		ions.density= density;
		Np_initial=(Nx-1)*dx*ions.density/ions.weight;	
		std::cout<<" Initial amount of i Particles simulated " << Np_initial <<std::endl;
		ions.Np=0;
		ions.Np_alloc=2*Np_initial;
		ions.part = new Particle[2*ions.Np_alloc];
		ions.n_avg = new double[Nx];
	
		//Initial conditions
		double v_thi = sqrt(2.0*Kb*ions.temperature*EV_TO_K/ions.mass);
	
		for (int i=0 ; i<Np_initial ;i++){
		
			x= sim.get_domain_start() + 0.5*rnd()*sim.get_domain_length();
			v= SampleVel2(v_thi);
			Add_Particle( &ions, x, v);
		}

		double v_the = sqrt(2.0*Kb*electrons.temperature*EV_TO_K/electrons.mass);
		for (int i=0 ; i<Np_initial ;i++){
		
			x= sim.get_domain_start() + 0.5*rnd()*sim.get_domain_length();
			v= SampleVel2(v_the);
			Add_Particle( &electrons, x, v);		
			}
	}
	else
	{
		load_Specie( &electrons,  0);
		load_Specie( &ions,  1);
		electrons.n_avg = new double[Nx];
		ions.n_avg = new double[Nx];

	}

	double v_thi = sqrt(2.0*Kb*ions.temperature*EV_TO_K/ions.mass);
	double v_the = sqrt(2.0*Kb*electrons.temperature*EV_TO_K/electrons.mass);

	sim.clean_charge();
	ScatterSpecie_2(&electrons, &sim);
	ScatterSpecie_2(&ions, &sim);

	if ( sim.solve_potential_CG(0.0, 0.0))
	{	}
	else std::cout<< " Potential not solved " <<std::endl;	
	
	sim.compute_ef();
	

	if(re_start == false)
	{
		RewindSpecies( &electrons, &sim, dt);
		RewindSpecies( &ions, &sim, dt);
	}		

	for (int i=0;i<Nx;i++){	 
		phi_avg[i]=0.0;
		electrons.n_avg[i]=0.0;
		ions.n_avg[i]=0.0;

	}


	int average_over=10000;

	for(int i_time=0;  i_time<=NUM_TS; i_time++)
	{	

		file_c<< i_time << "\t" << electrons.Np << "\t" << ions.Np <<std::endl;
		
		//Charge projection
		sim.clean_charge();
		ScatterSpecie_2(&electrons, &sim);
		ScatterSpecie_2(&ions, &sim);

		//Solve Poisson equation
		if ( !sim.solve_potential_CG(0.0, 0.0)) std::cout<< " Potential not solved " <<std::endl;	
		
		phi_aux= sim.get_potential();
		for(int ii=0; ii<Nx ;ii ++ ){
			phi_avg[ii] +=phi_aux[ii];
		}
		ScatterSpecies_den(&electrons, Nx ,dx, sim.get_domain_start());	
		ScatterSpecies_den(&ions, Nx ,dx, sim.get_domain_start());
		
		if(i_time%average_over==0){	
			for (int i=0;i<Nx;i++){	
				file_avg<< i*dx <<"\t"<< phi_avg[i]/(average_over*1.0)  <<"\t" << electrons.n_avg[i]/(average_over*1.0) << "\t" <<  ions.n_avg[i]/(average_over*1.0) <<std::endl; 
				phi_avg[i]=0.0;
				electrons.n_avg[i]=0.0;
				ions.n_avg[i]=0.0;

			}
			file_avg<< std::endl; 
			file_avg<< std::endl;
		}


		//Compute electric field
		sim.compute_ef();

		//Move particles
		Push_specie( &electrons, &sim, dt);
		Push_specie( &ions, &sim, dt);

		Reorder_Specie(&electrons);
		Reorder_Specie(&ions);
		thermalization( &electrons, &sim, v_the ,dt);

		reinjection_B(&electrons, &ions, &sim , v_the,  v_thi);		

		
		if(i_time%100000==0){	
			std::cout<<"iter = " << i_time << " Number of electrons = "<< electrons.Np << " Number of ions = "<< ions.Np  <<std::endl;			
			/*sim.print_charge_density();
			sim.print_potential();				
			sim.print_ef();	
			*/			
			save_Specie( &electrons,  0);
			save_Specie( &ions,  1);	
					
		}
		
		
	}	
	
	sim.~Domain();
	delete [] ions.part;
	delete [] electrons.part;			
	delete [] phi_aux;
	delete [] phi_avg;	
	delete [] electrons.n_avg;	
	delete [] ions.n_avg;
	file_c.close();
	
	return 0;


}
