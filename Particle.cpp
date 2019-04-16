#include <cmath>
#include "Particle.h"





Particle::Particle(){
	x=0.0;
	vx=0.0;
	alive=false;
}



void Particle::set_position(double xpos){
	x=xpos;
}

void Particle::set_velocity(double velx){
	vx=velx;
}
double Particle::get_position(void){
	return x;
}
double Particle::get_velocity(void){
	return vx;
}

bool Particle::is_alive(){
	return alive;
}

void Particle::set_is_alive(bool state){
	alive=state;
}

void Particle::move_particle_LF(double force, double dt){
	vx += force*dt;
	x += vx*dt;
}

void Particle::velocity_rewind_LF(double force, double dt){
	vx += force*dt;
}


Particle::~Particle(void){}
