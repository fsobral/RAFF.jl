/*

g++ -std=c++11 -I/usr/include/eigen3 -I/usr/local/include/theia/libraries/vlfeat -I/usr/include/suitesparse -I/usr/local/include/theia/libraries/ -I/usr/local/include/theia/libraries/statx -I/usr/local/include/theia/libraries/optimo -I/usr/local/include/theia ransac.cpp -o ransac -ltheia -lglog

 */

#include<theia.h>

// Our "data".
struct Point {
  double x; double y;
};

// Our "model".
struct Line {
  double m; double b;
};

// Estimator class.
class LineEstimator: public theia::Estimator<Point, Line> {

public:

   LineEstimator(){};

  ~LineEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override { return 2; }

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Line>* models) const override {
    Line model;
    model.m = (data[1].y - data[0].y)/(data[1].x - data[0].x);
    model.b = data[1].y - model.m*data[1].x;
    models->push_back(model);
    return true;
  }

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Line& line) const override {
    return point.y - (line.m*point.x + line.b);
  }
  
};

void run_ransac(void) {
  // Generate your input data using your desired method.
  // We put pseudo-code here for simplicity.
  std::vector<Point> input_data;

  // Add 700 inliers.
  for (int i = 0; i < 700; i++) {
    double x = i * 20.0 / 700.0;
    Point inlier_point = { x , 3.0 * x + 4.0 };
    input_data.push_back(inlier_point);
  }
  // Add 300 outliers.
  for (int i = 0; i < 300; i++) {
    double x = i * 20.0 / 300.0;
    Point outlier_point = { x , 4.0 };
    input_data.push_back(outlier_point);
  }

  // Specify RANSAC parameters.
  double error_threshold = 0.3;
  int min_num_inliers = 600;
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
