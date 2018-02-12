/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 200;
	default_random_engine gen;
	double std_x, std_y, std_theta;
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	weights = vector<double>(num_particles);
	particles = vector<Particle>(num_particles);
	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
		Particle particle;
		particle.id = i;
		particle.weight = 1.0;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		// not sure if I should append here the particle to the particles list
		weights[i] = 1.0; //not sure if I need this at this point
		is_initialized = true; // not sure about this line either
		particles[i] = particle;
	}
	cout << "initiation" << endl;
	cout << "gps_x: " << x << endl;
	cout << "gps_y: " << y << endl;
	cout << "gps_theta: " << theta << endl;
	cout << "particles[0].id: " << particles[0].id << endl;
	cout << "particles[0].x: " << particles[0].x << endl;
	cout << "particles[0].y: " << particles[0].y << endl;
	cout << "particles[0].theta: " << particles[0].theta << endl;
	cout << "particles[0].weight: " << particles[0].weight << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> x_noise(0, std_pos[0]);
	normal_distribution<double> y_noise(0, std_pos[1]);
	normal_distribution<double> theta_noise(0, std_pos[2]);
	default_random_engine gen;
	for (int i = 0; i < num_particles; i++) {
		Particle particle = particles[i];
		double theta = particle.theta;
		if (yaw_rate == 0) {
			particle.x += velocity * delta_t;
		}
		else {
			particle.x += (velocity / yaw_rate) * (sin(theta + yaw_rate * delta_t) - sin(theta)) + x_noise(gen);
			particle.y += (velocity / yaw_rate) * (cos(theta) - cos(theta + (yaw_rate * delta_t))) + y_noise(gen);
			particle.theta += yaw_rate * delta_t + theta_noise(gen);
		}
		// I'm not sure what noise I should add here, if at all
		particles[i] = particle;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < particles.size(); i++) {
		Particle particle = particles[i];
		particle.sense_x = std::vector<double>(observations.size());
		particle.sense_y = std::vector<double>(observations.size());
		particle.associations = std::vector<int>(observations.size());
		// TRANSFORMATION
		for (int j = 0; j < observations.size(); j++) {
			particle.sense_x[j] = particle.x + (cos(particle.theta) * observations[j].x) - (sin(particle.theta) * observations[j].y);
			particle.sense_y[j] = particle.y + (sin(particle.theta) * observations[j].x) + (cos(particle.theta) * observations[j].y);
			//particle.associations[j] = 0;
		}
		// ASSOCIATION
		for (int j = 0; j < particle.sense_x.size(); j++) {
			double prev_distance = 1000;
			particle.associations[j] = map_landmarks.landmark_list[0].id_i;
			for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
				Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
				double distance_x = particle.sense_x[j] - landmark.x_f;
				double distance_y = particle.sense_y[j] - landmark.y_f;
				double distance = sqrt(distance_x * distance_x + distance_y * distance_y);
				if (distance < prev_distance) {
					particle.associations[j] = landmark.id_i;
					prev_distance = distance;
				}
			}
		}
		double weight = 1.0;
		double gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));
		for (int j = 0; j < particle.associations.size(); j++) {			
			int landmark_id = particle.associations[j];
			LandmarkObs landmarkb = observations[j];
			Map::single_landmark_s landmark = map_landmarks.landmark_list[landmark_id - 1];
			double exponent = (pow(particle.sense_x[j] - landmark.x_f, 2)) / (2 * pow(std_landmark[0], 2)) + (pow(particle.sense_y[j] - landmark.y_f, 2)) / (2 * pow(std_landmark[1], 2));
			weight = weight * gauss_norm * exp(-exponent);
		}
		particle.weight = weight;
		particles[i] = particle;
	}
	double weights_sum = 0;
	for (int i = 0; i < particles.size(); i++) {
		weights_sum += particles[i].weight;
	}
	double efi = 0;
	for (int i = 0; i < particles.size(); i++) {
		particles[i].weight = particles[i].weight / weights_sum;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// NOTE: I'm trying a different kind of resampling method. hope it will work.

	vector<double> weights_agg(num_particles);
	double weights_sum = 0;
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < particles.size(); i++) {
		weights_sum = weights_sum + particles[i].weight;
		weights_agg[i] = weights_sum;
	}

	vector<Particle> new_particles(num_particles);
	for (int i = 0; i < particles.size(); i++) {
		// BINARY SEARCH
		double selected = distribution(generator);
		int top = num_particles;
		int bottom = 0;
		int range = top - bottom;
		int to_check = top - (range / 2);
		double cell = weights_agg[to_check];
		for (int i = 0; i < ceil(log2(num_particles)); i++) {
			if (cell < selected) {
				bottom = to_check;
				range = top - bottom;
				to_check = top - (range / 2);
				cell = weights_agg[to_check];
				if (range == 1) {
					break;
				}
			}
			if (cell > selected) {
				top = to_check;
				range = top - bottom;
				to_check = top - (range / 2);
				cell = weights_agg[to_check];
				if (range == 1) {
					break;
				}
			}
			if (cell == selected) {
				break;
			}
		}
		new_particles[i] = particles[to_check];
		//double x = new_particles[i].x;
		//double y = new_particles[i].y;
		//double theta = new_particles[i].theta;
		//default_random_engine gen;
		//double std_x, std_y, std_theta;
		//std_x = sigma_landmark[0];
		//std_y = sigma_landmark[1];
		//std_theta = 0.1;
		//normal_distribution<double> dist_x(x, std_x);
		//normal_distribution<double> dist_y(y, std_y);
		//normal_distribution<double> dist_theta(theta, std_theta);
		//new_particles[i].x = dist_x(gen);
		//new_particles[i].y = dist_y(gen);
		//new_particles[i].theta = dist_theta(gen);
	}
	particles = new_particles;
	weights_sum = 0;
	for (int i = 0; i < particles.size(); i++) {
		weights_sum += particles[i].weight;
	}
	int counter = 0;
	for (int i = 0; i < particles.size(); i++) {
		particles[i].weight = particles[i].weight / weights_sum;
		if (particles[i].weight > 0.01) {
			counter += 1;
		}
	}
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
