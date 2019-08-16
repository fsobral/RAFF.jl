/*

g++ -std=c++11 -I/usr/include/eigen3 -I/usr/local/include/theia/libraries/vlfeat -I/usr/include/suitesparse -I/usr/local/include/theia/libraries/ -I/usr/local/include/theia/libraries/statx -I/usr/local/include/theia/libraries/optimo -I/usr/local/include/theia ransac.cpp -o ransac -ltheia -lglog -lceres

 */

#include "ransac.hpp"
#include "ceres/ceres.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::SoftLOneLoss;

/* Implementation of the line estimator. */

double LineEstimator::SampleSize() const { return 5; }

bool LineEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Line>* models) const {

  Problem problem;

  Line model;

  double a = -100.0;
  double b = 2000.0;

  // Use ceres to solve LS problems
  
  for (int i = 0; i < this->SampleSize(); ++i) {

    problem.AddResidualBlock(
        new AutoDiffCostFunction<LinearResidual, 1, 1, 1>(
            new LinearResidual(data[i].x, data[i].y)),
        new SoftLOneLoss(0.5),
        &a, &b);

  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Add the results to Theia

  model.a = a;
  model.b = b;

  models->push_back(model);
  
  return true;
}

double LineEstimator::Error(const Point& point, const Line& line) const {
  return point.y - (line.a * point.x + line.b);
}

//template <typename T>
void run_ransac(int argc, char** argv) {
  // Generate your input data using your desired method.
  // We put pseudo-code here for simplicity.
  std::vector<Point> input_data;

  double x, y, z;
  
  // Add data
  std::ifstream file("/tmp/output.txt");

  // Discard first line of file
  file >> x;
  
  while (file >> x >> y >> z) {

    Point point = {x, y};
    
    input_data.push_back(point);
    
    //std::cout << x << ' ' << y << std::endl;

  }

  file.close();
  
  // Estimate the line with RANSAC.
  LineEstimator line_estimator;
  Line best_line;

  // Set the ransac parameters.
  theia::RansacParameters params;
  params.error_thresh = 0.1;
  params.use_mle = true;
  params.min_inlier_ratio = 0.5;

  // Create Ransac object, specifying the number of points to sample to
  // generate a model estimation.
  theia::Ransac<LineEstimator> ransac_estimator(params, line_estimator);
  // Initialize must always be called!
  ransac_estimator.Initialize();

  theia::RansacSummary summary;
  ransac_estimator.Estimate(input_data, &best_line, &summary);

  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(15);
  std::cout << "RANSAC:" << std::setw(25) << best_line.a << "," << std::setw(25) << best_line.b << std::endl;

  // Create Prosac object
  theia::Prosac<LineEstimator> prosac_estimator(params, line_estimator);
  // Initialize must always be called!
  prosac_estimator.Initialize();

  prosac_estimator.Estimate(input_data, &best_line, &summary);

  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(15);
  std::cout << "PROSAC:" << std::setw(25) << best_line.a << "," << std::setw(25) << best_line.b << std::endl;
  
}

int main(int argc, char** argv) {

  run_ransac(argc, argv);

  return 0;

}
