/*

g++ -std=c++11 -I/usr/include/eigen3 -I/usr/local/include/theia/libraries/vlfeat -I/usr/include/suitesparse -I/usr/local/include/theia/libraries/ -I/usr/local/include/theia/libraries/statx -I/usr/local/include/theia/libraries/optimo -I/usr/local/include/theia ransac.cpp -o ransac -ltheia -lglog

 */

#include "ransac.hpp"
#include <fstream>
#include <iostream>

double LineEstimator::SampleSize() const { return 2; }

bool LineEstimator::EstimateModel(const std::vector<Point>& data,
                     std::vector<Line>* models) const {
  Line model;
  model.m = (data[1].y - data[0].y)/(data[1].x - data[0].x);
  model.b = data[1].y - model.m*data[1].x;
  models->push_back(model);
  return true;
}

double LineEstimator::Error(const Point& point, const Line& line) const {
  return point.y - (line.m*point.x + line.b);
}

void run_ransac(void) {
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
    
    std::cout << x << ' ' << y << std::endl;

  }

  file.close();
  
  // Specify RANSAC parameters.
  double error_threshold = 0.3;
  int min_num_inliers = 10;
  int max_iters = 1000;

  // Estimate the line with RANSAC.
  LineEstimator line_estimator;
  Line best_line;
  // Set the ransac parameters.
  theia::RansacParameters params;
  params.error_thresh = 0.1;

  // Create Ransac object, specifying the number of points to sample to
  // generate a model estimation.
  theia::Ransac<LineEstimator> ransac_estimator(params, line_estimator);
  // Initialize must always be called!
  ransac_estimator.Initialize();

  theia::RansacSummary summary;
  ransac_estimator.Estimate(input_data, &best_line, &summary);
  LOG(INFO) << "Line m = " << best_line.m << "*x + " << best_line.b;

}

int main(void) {

  run_ransac();

  return 0;

}
