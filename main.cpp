#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <iostream>
#include <assert.h>

namespace sirAlgorithm {

const double xInitialValue_ = 0.1;
const double varianceY_ = 10;
const double varianceX_ = 1;
size_t numberParticles_;
size_t timestep_ = 1;
std::default_random_engine generator_;
std::normal_distribution<double> xNoise_(0.0, varianceX_);
std::normal_distribution<double> yNoise_(0.0, varianceY_);
std::vector<double> particleWeights_;

double simulateXParticle(double previousXValue) {
  double newParticle = 0.5*previousXValue+25*previousXValue/(1+std::pow(previousXValue, 2))+
      8*std::cos(1.2*timestep_) + xNoise_(generator_);
  return newParticle;
}

double simulateYParticle(double currentXValue) {
  return std::pow(currentXValue, 2)/20 + yNoise_(generator_);
}

std::vector<double> updateX_k(std::vector<double> xPreviousValues){
  std::vector<double> newXValues;
  if (!xPreviousValues.size()) {
    for (size_t particle = 0; particle < numberParticles_; particle++) {
      newXValues.push_back(simulateXParticle(xInitialValue_));
    }
  } else {
    assert(xPreviousValues.size() == numberParticles_);
    for (size_t particle = 0; particle < numberParticles_; particle++) {
      newXValues.push_back(simulateXParticle(xPreviousValues.at(particle)));
    }
  }
  return newXValues;
}
std::vector<double> updateY_k(std::vector<double> xCurrentValues){
  std::vector<double> newYValues;
  for (size_t particle = 0; particle < numberParticles_; particle++) {
    newYValues.push_back(simulateYParticle(xCurrentValues.at(particle)));
  }
  return newYValues;
}
std::vector<double> updateWeights(){}
double calculateYLikelihood(){}
std::vector<double> resampleParticles(){}

int main()
{
  std::vector<double> x_k;
  std::vector<double> y_k;
  std::cout << "Number of particles: ";
  std::cin >> numberParticles_;
  assert(numberParticles_ > 0);
  double particleInitialWeights = 1.0/numberParticles_;
  particleWeights_.assign(numberParticles_, particleInitialWeights);
  std::cout << "\nNumber of time steps: ";
  int timeStepCount;
  std::cin >> timeStepCount;
  assert(timeStepCount > 0);

  return 0;
}
} // namespace sirAlgorithm
