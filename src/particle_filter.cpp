/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
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

  // number of particles
  num_particles = 100;

  // creating normal (Gaussian) distributions for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // initial weight
  const double new_weight = 1.0;

  // let's initialize all of them, then
  for ( int i = 0;  i < num_particles;  ++i ) {

    // creating a new particle
    Particle new_particle = {
      -1,                      // "id", khm
      dist_x(generator),      // sampled x value
      dist_y(generator),      // sampled y value
      dist_theta(generator),  // sampled theta value
      new_weight              // initial weight of the particle
    };

    // filling both the "weights" and the "particles" vectors 
    weights.push_back(new_weight);
    particles.push_back(new_particle);
  }

  // we're all done
  is_initialized = true;

  // ...but the initialization of our internal landmark mapping lies still ahead
  is_landmark_mapping_initialized = false;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // defining constant terms outside the for loop
  const double epsilon = 0.001;
  const double c = abs(yaw_rate) < epsilon ? 1.0 : velocity / yaw_rate;
  const double distance = velocity * delta_t;
  const double delta_theta = yaw_rate * delta_t;

  // processing all the particles
  for ( int i = 0;  i < num_particles;  ++i ) {

    // implementing two different motion models depending on the yaw rate value
    if ( abs(yaw_rate) < epsilon ) {

      // using the basic motion model (theta remains unchanged)
      particles[i].x += cos(particles[i].theta) * distance;
      particles[i].y += sin(particles[i].theta) * distance;

    } else {

      // using the slightly more complicated motion model for a turning car 
      double new_theta = particles[i].theta + delta_theta;

      particles[i].x += c * (sin(new_theta) - sin(particles[i].theta));
      particles[i].y += c * (cos(particles[i].theta) - cos(new_theta));
      particles[i].theta = new_theta;
    }

    // creating normal (Gaussian) distributions accounting for measurement noise
    // around the calculated values acting as means
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

    // setting particle position -- including measurement noise
    particles[i].x = dist_x(generator);
    particles[i].y = dist_y(generator);
    particles[i].theta = dist_theta(generator);
  }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {

	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  // ?.... I'm not sure whether 'predicted' expresses the right concept here.
  // Whatever. I pass landmarks with fixed positions to the function, anyway.

  // defining the variables only once
  double distance, best_distance;
  size_t i, j;
  const size_t observation_count = observations.size();
  const size_t landmark_count = predicted.size();

  for( i = 0;  i < observation_count;  ++i ) {

    // initializing helper variable for this observation
    best_distance = 9999.0;

    for ( j = 0;  j < landmark_count;  ++j ) {

      // measuring the distance between the landmark and the observation
      distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      // in case this is a better value compared to what we've already seen so far,
      // we'll mark it as the best one
      if ( distance < best_distance ) {

        best_distance = distance;
        observations[i].id = predicted[j].id;
      }
    }
  }
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

  // helper variables
  double theta, landmark_dist, weight, p;
  std::vector<LandmarkObs> transformed_observations;  // observations in map coordinates
  std::vector<LandmarkObs> landmarks_in_range;        // available landmarks within sensor range for the current particle
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;

  // more helpers (avoiding function call in loop condition checking -- maybe the compiler would do it anyway, too)
  const size_t observation_count = observations.size();
  const size_t map_landmark_count = map_landmarks.landmark_list.size();
  size_t j;

  LandmarkObs observation;                            // observation in map coordinates
  LandmarkObs landmark;                               // a landmark within a given particle's sensor range
  LandmarkObs nearest_landmark;                       // nearest landmark around an observation

  // initializing the associative array for the landmarks (doing this only once)
  if ( !is_landmark_mapping_initialized ) {
    initLandmarkMapping(map_landmarks.landmark_list);
  }

  // processing all the particles
  for ( int i = 0;  i < num_particles;  ++i ) {

    // this one is only for enhanced readability
    theta = particles[i].theta;

    // initializing helpers 
    transformed_observations.clear();
    landmarks_in_range.clear();
    associations.clear();
    sense_x.clear();
    sense_y.clear();

    // transforming observation to map coordinates
    for ( j = 0;  j < observation_count;  ++j ) {

      // performing the homogenous transformation for the particle
      observation.x = particles[i].x + cos(theta) * observations[j].x - sin(theta) * observations[j].y;
      observation.y = particles[i].y + sin(theta) * observations[j].x + cos(theta) * observations[j].y;

      transformed_observations.push_back(observation);
    }

    // collecting available landmarks around the particle within the sensor range
    for( j = 0;  j < map_landmark_count;  ++j ) {

      // storing a landmark in case it's within the current particle's range
      landmark_dist = dist(particles[i].x, particles[i].y, (double)map_landmarks.landmark_list[j].x_f, (double)map_landmarks.landmark_list[j].y_f);
      if ( landmark_dist <= sensor_range ) {

        landmark = {
          map_landmarks.landmark_list[j].id_i,
          (double)map_landmarks.landmark_list[j].x_f,
          (double)map_landmarks.landmark_list[j].y_f
        };

        // storing important values
        landmarks_in_range.push_back(landmark);
        associations.push_back(landmark.id);
        sense_x.push_back(landmark.x);
        sense_y.push_back(landmark.y);
      }
    }

    // setting associations (...)
    SetAssociations(particles[i], associations, sense_x, sense_y);

    // performing nearest neighbor association
    dataAssociation(landmarks_in_range, transformed_observations);

    // this variable will hold the new weight of the particle
    weight = 1.0;

    // processing the result
    for ( j = 0;  j < observation_count;  ++j ) {

      // getting the position of the chosen landmark in a simple way
      // (that's why we prepared the member variable 'landmark_mapping')
      nearest_landmark = landmark_mapping[transformed_observations[j].id];

      // we've found the chosen landmark (the following function is defined in helper_functions.h)
      p = multivariateGaussian(transformed_observations[j].x, transformed_observations[j].y,
                               nearest_landmark.x, nearest_landmark.y,
                               std_landmark[0], std_landmark[1]);
      weight *= p;
    }

    // storing the new weight
    weights[i] = weight;
    particles[i].weight = weight;
  }
}


/**
 * initLandmarkMapping -- preparing an associative array (=a C++ map) with the landmark IDs as keys
 * and LandmarkObs as values so that the coordinates of a given landmark can be retrieved in O(1).
 * @param map_landmark_list: a vector of 'single_landmark_s' data
 */
void ParticleFilter::initLandmarkMapping(const std::vector<Map::single_landmark_s> &map_landmark_list) {

  // trying to optimize code efficiency (I'm not sure whether I really need this, though :-/ )  
  const size_t map_landmark_count = map_landmark_list.size();
  LandmarkObs landmark;

  // looping through all landmarks and making them available through a map with the IDs as keys
  for( size_t i = 0;  i < map_landmark_count;  ++i ) {

    landmark = {
      map_landmark_list[i].id_i,
      (double)map_landmark_list[i].x_f,
      (double)map_landmark_list[i].y_f
    };

    landmark_mapping[map_landmark_list[i].id_i] = landmark;
  }

  is_landmark_mapping_initialized = true; 
}



void ParticleFilter::resample() {

	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // these variables will hold the new particles and weights data
    std::vector<Particle> resampled_particles;
    std::vector<double> resampled_weights;

    // helper variables
    double beta = 0.0;
    auto max_weight_iterator = std::max_element(weights.begin(), weights.end());
    double max_weight = *max_weight_iterator;

    // getting random numbers from a uniform distibution
    std::uniform_real_distribution<double> dist_uniform(0.0, 2.0 * max_weight);

    // ...and getting the first random weight index
    std::discrete_distribution<unsigned int> dist_index(weights.begin(), weights.end());
    unsigned int index = dist_index(generator);

    // processing the particles
    for ( int i = 0;  i < num_particles; ++i ) {

       beta += dist_uniform(generator);

       // looking for the next candidate
       while(weights[index] < beta) {
         beta -= weights[index];
         index = (index + 1) % num_particles;
       }

      resampled_particles.push_back(particles[index]);
      resampled_weights.push_back(weights[index]);
    }

    // setting the new world order
    particles = resampled_particles;
    weights = resampled_weights;
}


void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    // particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations = associations;
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
