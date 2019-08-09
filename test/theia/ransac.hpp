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
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Line>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Line& line) const override;

};
  
void run_ransac(void);
