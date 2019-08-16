#include<theia.h>

// Our "data".
struct Point {
  double x; double y;
};

// Our "model".
struct Line {
  double a; double b;
};

// Structure to model LS solver for Ceres. This structure is to find
// linear one-dimensinal models.
struct LinearResidual {
  LinearResidual(double x, double y)
      : x_(x), y_(y) {}

  template <typename T> bool operator()(const T* const a,
                                        const T* const b,
                                        T* residual) const {
    residual[0] = y_ - (a[0] * x_ + b[0]);
    return true;
  }

 private:
  const double x_;
  const double y_;
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
