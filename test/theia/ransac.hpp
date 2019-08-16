#include<theia.h>
#include<math.h>

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

// Exponential model

// Our "model".
struct Exponential {
  double a; double b; double c;
};

// Structure to model LS solver for Ceres. This structure is to find
// linear one-dimensinal models.
struct ExponentialResidual {
  ExponentialResidual(double x, double y)
      : x_(x), y_(y) {}

  template <typename T> bool operator()(const T* const p,
                                        T* residual) const {
    residual[0] = y_ - (p[0] + p[1] * exp(- p[2] * x_));
    return true;
  }

 private:
  const double x_;
  const double y_;
};

// Estimator class.
class ExponentialEstimator: public theia::Estimator<Point, Exponential> {

public:

   ExponentialEstimator(){};

  ~ExponentialEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Exponential>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Exponential& line) const override;

};

void run_ransac(void);
