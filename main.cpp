#define _USE_MATH_DEFINES

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include <assert.h>
#include <math.h>

constexpr double xInitialValue_ = 0.1;
constexpr double varianceY_ = 10;
constexpr double varianceX_ = 1;
size_t numberParticles_;
size_t timestep_ = 1;
std::default_random_engine generator_;
std::normal_distribution<double> xNoise_(0.0, varianceX_);
std::normal_distribution<double> yNoise_(0.0, varianceY_);

double simulateXParticle(double previousXValue) {
  double newParticle = 0.5*previousXValue+25*previousXValue/(1+std::pow(previousXValue, 2))+
      8*std::cos(1.2*timestep_) + xNoise_(generator_);
  return newParticle;
}

double simulateYParticle(double currentXValue, double* yParticleNoise) {
  *yParticleNoise = yNoise_(generator_);
  return std::pow(currentXValue, 2)/20 + *yParticleNoise;
}

void updateX_k(std::vector<double> &xPreviousValues){
  assert(xPreviousValues.size() == numberParticles_);
  for (size_t particle = 0; particle < numberParticles_; particle++) {
    xPreviousValues.at(particle) = simulateXParticle(xPreviousValues.at(particle));
  }
}
std::vector<double> updateY_k(std::vector<double> &xCurrentValues, double* yParticleNoise){
  std::vector<double> newYValues;
  for (size_t particle = 0; particle < numberParticles_; particle++) {
    newYValues.push_back(simulateYParticle(xCurrentValues.at(particle), yParticleNoise);
    double yParticleLikelihood = estimateLikelihood(yParticleNoise);
    updateWeights();
  }
  return newYValues;
}

double estimateLikelihood(double* yParticleNoise) {
  double likelihood = std::sqrt(0.5*M_1_PI*1/varianceY_);
  likelihood *= std::exp(-1/(2*varianceY_)*std::pow(*yParticleNoise, 2));
  return likelihood;
}

void updateWeights(std::vector<double> &weightsVector,
                   const std::vector<double> &yVector,
                   const std::vector<double> &xVector){
  for (size_t particle = 0; particle < numberParticles_; particle++) {
    double likelihood = estimateLikelihood(yVector.at(particle), xVector.at(particle));
    weightsVector.at(particle) *= likelihood;
  }
  double weightSum = std::accumulate(weightsVector.begin(), weightsVector.end(), 0);
  auto normalize = [&weightSum](const double& n) { return n/weightSum; };
  std::for_each(weightsVector.begin(), weightsVector.end(), normalize);
}
double calculateYLikelihood(){}
std::vector<double> resampleParticles(){}

int main()
{
  std::unique_ptr<double> yParticleNoise(new double);
  std::vector<double> y_k;
  std::cout << "Number of particles: ";
  std::cin >> numberParticles_;
  assert(numberParticles_ > 0);
  double particleInitialWeights = 1.0/numberParticles_;
  std::vector<double> particleWeights_(numberParticles_, particleInitialWeights);
  std::vector<double> x_k(numberParticles_, xInitialValue_);
  updateX_k(x_k);
  updateYandWeights(x_k, yParticleNoise);
  std::cout << "\nNumber of time steps: ";
  int timeStepCount;
  std::cin >> timeStepCount;
  assert(timeStepCount > 0);

  return 0;
}
