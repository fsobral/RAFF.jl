/*

g++ -std=c++11 -I/usr/include/eigen3 -I/usr/local/include/theia/libraries/vlfeat -I/usr/include/suitesparse -I/usr/local/include/theia/libraries/ -I/usr/local/include/theia/libraries/statx -I/usr/local/include/theia/libraries/optimo -I/usr/local/include/theia -L/usr/local/lib ransac.cpp -o ransac -ltheia -lglog -lceres

 */

#include "ransac.hpp"
#include "ceres/ceres.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::SoftLOneLoss;

/* ************************************* */
/* Implementation of the line estimator. */
/* ************************************* */

std::ostream & operator<<(std::ostream & Str, Line const & v) { 
  Str << std::setw(25) << v.a << " " << std::setw(25) << v.b;
  return Str;
}

LineEstimator::LineEstimator(int ssize_) {ssize = ssize_;}

double LineEstimator::SampleSize() const { return ssize; }

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
            new LinearResidual(data[i].x[0], data[i].y)),
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
  return point.y - (line.a * point.x[0] + line.b);
}

/* ************************************** */
/* Implementation of the cubic estimator. */
/* ************************************** */

std::ostream & operator<<(std::ostream & Str, Cubic const & v) { 
  Str << std::setw(25) << v.a << " " << std::setw(25) << v.b << " " <<
    std::setw(25) << v.c << " " << std::setw(25) << v.d;
  return Str;
}

CubicEstimator::CubicEstimator(int ssize_) {ssize = ssize_;}

double CubicEstimator::SampleSize() const { return ssize; }

bool CubicEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Cubic>* models) const {

  Problem problem;

  Cubic model;

  double p[] = {0.0, 0.0, 0.0, 0.0};

  // Use ceres to solve LS problems
  
  for (int i = 0; i < this->SampleSize(); ++i) {

    problem.AddResidualBlock(
        new AutoDiffCostFunction<CubicResidual, 1, 4>(
            new CubicResidual(data[i].x[0], data[i].y)),
        new SoftLOneLoss(0.5),
        p);

  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Add the results to Theia

  model.a = p[0];
  model.b = p[1];
  model.c = p[2];
  model.d = p[3];

  models->push_back(model);
  
  return true;
}

double CubicEstimator::Error(const Point& point, const Cubic& p) const {
  double px = point.x[0];
  return point.y - (p.a * std::pow(px, 3) + p.b * std::pow(px, 2) + p.c * px + p.d);
}

/* ******************************************** */
/* Implementation of the exponential estimator. */
/* ******************************************** */

std::ostream & operator<<(std::ostream & Str, Exponential const & v) { 
  Str << std::setw(25) << v.a << " " << std::setw(25) << v.b << " " << std::setw(25) << v.c;
  return Str;
}

ExponentialEstimator::ExponentialEstimator(int ssize_) {ssize = ssize_;}

double ExponentialEstimator::SampleSize() const { return ssize; }

bool ExponentialEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Exponential>* models) const {

  Problem problem;

  Exponential model;

  double p[] = {0.0, 0.0, 0.0};

  // Use ceres to solve LS problems
  
  for (int i = 0; i < this->SampleSize(); ++i) {

    problem.AddResidualBlock(
        new AutoDiffCostFunction<ExponentialResidual, 1, 3>(
            new ExponentialResidual(data[i].x[0], data[i].y)),
        new SoftLOneLoss(0.5),
        p);

  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Add the results to Theia

  model.a = p[0];
  model.b = p[1];
  model.c = p[2];

  models->push_back(model);
  
  return true;
}

double ExponentialEstimator::Error(const Point& point, const Exponential& p) const {
  return point.y - (p.a + p.b * std::exp(- p.c * point.x[0]));
}

/* ***************************************** */
/* Implementation of the logistic estimator. */
/* ***************************************** */

std::ostream & operator<<(std::ostream & Str, Logistic const & v) { 
  Str << std::setw(25) << v.a << " " << std::setw(25) << v.b << " " <<
    std::setw(25) << v.c << " " << std::setw(25) << v.d;
  return Str;
}

LogisticEstimator::LogisticEstimator(int ssize_) {ssize = ssize_;}

double LogisticEstimator::SampleSize() const { return ssize; }

bool LogisticEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Logistic>* models) const {

  Problem problem;

  Logistic model;

  double p[] = {0.0, 0.0, 0.0, 0.0};

  // Use ceres to solve LS problems
  
  for (int i = 0; i < this->SampleSize(); ++i) {

    problem.AddResidualBlock(
        new AutoDiffCostFunction<LogisticResidual, 1, 4>(
            new LogisticResidual(data[i].x[0], data[i].y)),
        new SoftLOneLoss(0.5),
        p);

  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Add the results to Theia

  model.a = p[0];
  model.b = p[1];
  model.c = p[2];
  model.d = p[3];

  models->push_back(model);
  
  return true;
}

double LogisticEstimator::Error(const Point& point, const Logistic& p) const {
  return point.y - (p.a + p.b / (1.0 + std::exp(- p.c * point.x[0] + p.d)));
}

/* *************************************** */
/* Implementation of the circle estimator. */
/* *************************************** */

std::ostream & operator<<(std::ostream & Str, Circle const & v) { 
  Str << std::setw(25) << v.x << " " << std::setw(25) << v.y << " " << std::setw(25) << v.r;
  return Str;
}

double CircleEstimator::SampleSize() const { return ssize; }

bool CircleEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Circle>* models) const {

  Problem problem;

  Circle model;

  double p[] = {0.0, 0.0, 1.0};

  // Use ceres to solve LS problems
  
  for (int i = 0; i < this->SampleSize(); ++i) {

    problem.AddResidualBlock(
        new AutoDiffCostFunction<CircleResidual, 1, 3>(
          new CircleResidual(data[i].x[0], data[i].x[1], data[i].y)),
        new SoftLOneLoss(0.5),
        p);

  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Add the results to Theia

  //std::cout << p[1] << ' ' << p[2] << std::endl;
  
  model.x = p[0];
  model.y = p[1];
  model.r = p[2];

  models->push_back(model);
  
  return true;
}

CircleEstimator::CircleEstimator(int ssize_) {ssize = ssize_;}

double CircleEstimator::Error(const Point& point, const Circle& p) const {
  return point.y - (std::pow(point.x[0] - p.x, 2) + std::pow(point.x[1] - p.y, 2) - p.r * p.r);
}

/*
  Compute the quality measure.
*/

template <class T, class M>
double ls_measure(T estimator, M solution, std::vector<Point> input_data,
		  std::vector<bool> outliers){

  double sr = 0.0;
  double err;

  for (int i = 0; i < input_data.size(); i++) {

    if (!outliers[i]) {

      err = estimator.Error(input_data[i], solution);

      sr = sr + err * err;

    }

  }

  return sr;
  
}

/* Caller */

template <class T, class M>
void run_ransac(std::string filename, double error_thres, double min_inlier_ratio,
		double use_mle, double ftrusted, int max_iter, int multistart) {

  std::vector<Point> input_data;
  std::vector<bool> outliers;

  double y, z;
  int N;

  std::clock_t tini, tend;
  
  // Add data
  std::ifstream file(filename);

  // Get dimension of x
  file >> N;
  
  while (file) {

    std::vector<double> x(N, 0.0);

    if (N == 1) {

      file >> x[0] >> y >> z;
      if (x[0] == 0.0 && y == 0.0 && z == 0.0) break;

    }
    else if (N == 2) {

      file >> x[0] >> x[1] >> y >> z;
      if (x[0] == 0.0 && x[1] == 0.0 && y == 0.0 && z == 0.0) break;

    }
      
    Point point = {x, y};
    
    input_data.push_back(point);

    outliers.push_back((z != 0.0));
    
    // std::cout << x[0] << ' ' << y << std::endl;

  }

  file.close();

  T estimator(std::max(2, (int) (ftrusted * input_data.size())));
  M best_model;

  // Set the ransac parameters.
  theia::RansacParameters params;
  params.error_thresh = error_thres;
  params.use_mle = use_mle;
  params.min_inlier_ratio = min_inlier_ratio;
  params.max_iterations = max_iter;

  theia::RansacSummary summary;

  // Create Ransac object, specifying the number of points to sample to
  // generate a model estimation.

  for (int i = 0; i < multistart; i++) {
  
    tini = std::clock();
  
    theia::Ransac<T> ransac_estimator(params, estimator);
    // Initialize must always be called!
    ransac_estimator.Initialize();

    ransac_estimator.Estimate(input_data, &best_model, &summary);

    tend = std::clock();
  
    std::cout.unsetf(std::ios::scientific);
    std::cout.setf(std::ios::fixed);
    std::cout << std::setw(10) << "RANSAC";
    std::cout << std::setw(10) << std::setprecision(5) << (tend - tini) / (1.0 * CLOCKS_PER_SEC) << " ";
    std::cout << std::setw(10) << summary.num_iterations << " ";
    std::cout.unsetf(std::ios::fixed);
    std::cout.setf(std::ios::scientific);
    std::cout << std::setw(20) << std::setprecision(15);
    std::cout << ls_measure(estimator, best_model, input_data, outliers) << " ";
    std::cout << best_model << std::endl;

  }
  
  // // Create Prosac object
  // // This estimator assumes that the data is sorted by quality.

  // tini = std::clock();
  
  // theia::Prosac<T> prosac_estimator(params, estimator);
  // // Initialize must always be called!
  // prosac_estimator.Initialize();

  // prosac_estimator.Estimate(input_data, &best_model, &summary);

  // tend = std::clock();

  // std::cout.unsetf(std::ios::fixed);
  // std::cout.setf(std::ios::scientific);
  // std::cout << std::setprecision(15);
  // std::cout << std::setw(10) << "PROSAC:" << best_model << " & ";
  // std::cout.unsetf(std::ios::scientific);
  // std::cout.setf(std::ios::fixed);
  // std::cout << std::setw(10) << std::setprecision(5) << (tend - tini) / (1.0 * CLOCKS_PER_SEC) << std::endl;
  
  // Create LMed object

  for (int i = 0; i < multistart; i++) {

    tini = std::clock();
  
    theia::LMed<T> lmed_estimator(params, estimator);
    // Initialize must always be called!
    lmed_estimator.Initialize();

    lmed_estimator.Estimate(input_data, &best_model, &summary);

    tend = std::clock();

    std::cout.unsetf(std::ios::scientific);
    std::cout.setf(std::ios::fixed);
    std::cout << std::setw(10) << "LMED";
    std::cout << std::setw(10) << std::setprecision(5) << (tend - tini) / (1.0 * CLOCKS_PER_SEC) << " ";
    std::cout << std::setw(10) << summary.num_iterations << " ";
    std::cout.unsetf(std::ios::fixed);
    std::cout.setf(std::ios::scientific);
    std::cout << std::setw(20) << std::setprecision(15);
    std::cout << ls_measure(estimator, best_model, input_data, outliers) << " ";
    std::cout << best_model << std::endl;

  }
  
}

int main(int argc, char** argv) {

  double error_thres;
  double min_inlier_ratio = 0.0;
  bool use_mle = true;
  double ftrusted = 0.5;
  int max_iter = 1000;
  int ms = 1;
  std::string filename = std::string("/tmp/output.txt");

  if (argc < 3) {

    std::cout << "Call: " << argv[0] << " model err_thres [-mir mir] [-nomle] [-ft ft] [-f filename] [-max max_iter] [-ms multstarts]" << std::endl;

    return 1;

  }

  std::string s = std::string(argv[1]);
  error_thres = std::stod(argv[2]);
  
  int i = 3;

  while (i < argc)  {

    std::string tmp = std::string(argv[i]);

    if (tmp == "-mir") min_inlier_ratio = std::stod(argv[++i]);
    else if (tmp == "-nomle") use_mle = false;
    else if (tmp == "-ft") ftrusted = std::stod(argv[++i]);
    else if (tmp == "-f") filename = std::string(argv[++i]);
    else if (tmp == "-max") max_iter = std::stoi(argv[++i]);
    else if (tmp == "-ms") ms = std::stoi(argv[++i]);
    else {

      std::cout << "Call: " << argv[0] << " err_thres [-mir mir] [-mle]" << std::endl;
      return 1;

    }

    i += 1;

  }
  
  if (s == "linear") run_ransac<LineEstimator, Line>(filename, error_thres, min_inlier_ratio, use_mle, ftrusted, max_iter, ms);
  else if (s == "cubic") run_ransac<CubicEstimator, Cubic>(filename, error_thres, min_inlier_ratio, use_mle, ftrusted, max_iter, ms);
  else if (s == "expon") run_ransac<ExponentialEstimator, Exponential>(filename, error_thres, min_inlier_ratio, use_mle, ftrusted, max_iter, ms);
  else if (s == "logistic") run_ransac<LogisticEstimator, Logistic>(filename, error_thres, min_inlier_ratio, use_mle, ftrusted, max_iter, ms);
  else if (s == "circle") run_ransac<CircleEstimator, Circle>(filename, error_thres, min_inlier_ratio, use_mle, ftrusted, max_iter, ms);
  else std::cout << "Unknown model." << std::endl;

  return 0;

}
