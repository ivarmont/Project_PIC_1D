
#ifndef __PARTICLE_CLASS_INCLUDED__
#define __PARTICLE_CLASS_INCLUDED__
class Particle{

	private:
		double x,vx;
		bool alive;
		double mass, charge, weight;

	public:
		Particle();	
		void set_position(double);
		void set_velocity(double);
		double get_position(void);
		double get_velocity(void);
		void set_is_alive(bool);
		bool is_alive();
		void move_particle_LF(double, double);
		~Particle();
	
};

#endif
