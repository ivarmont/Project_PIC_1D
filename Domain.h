
#ifndef __DOMAIN_CLASS__
#define __DOMAIN_CLASS__

class Domain{
	private:
		int Nx;
		double dx;
		double x0;
		double xl;
		double x_s0;	
		double x_s1;
		/*
		double *n_e;
		double *n_pi;		
		double *n_ni;	
		*/
		double *rho;
		double *phi;
		double *ef;

		//for CG solver
		double *res;
		double *cg_dir;
		double *Acg_dir;
		

	public:
		Domain();
		Domain(int, double);	
		~Domain();
		double get_domain_length(void);
		double get_domain_start(void);
		double get_domain_end(void);
		
		double get_domain_source_begin(void);		
		double get_domain_source_end(void);	
		
		void set_charge(double, double);
		void set_charge_2(double, double);
		
		void clean_charge(void);
		
		bool solve_potential_GS(double, double);
		bool solve_potential_CG(double, double);		
		
		void compute_ef();

		double interpolate_ef(double);
		double interpolate_ef_2(double);

		double * get_potential();

		void print_charge_density();
		void print_potential();
		void print_ef();
};


#endif
