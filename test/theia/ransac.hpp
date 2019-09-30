#include<theia.h>
#include<math.h>

// Our "data".
struct Point {
  std::vector<double> x; double y;
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

  LineEstimator(int ssize_);

  ~LineEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Line>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Line& line) const override;

  private:

  int ssize;

};

/* ************ */
/* Cubic model. */
/* ************ */

// Our "model".
struct Cubic {
  double a; double b; double c; double d;
};

// Structure to model LS solver for Ceres. This structure is to find
// linear one-dimensinal models.
struct CubicResidual {
  CubicResidual(double x, double y)
      : x_(x), y_(y) {}

  template <typename T> bool operator()(const T* const p,
                                        T* residual) const {
    residual[0] = y_ - (((p[0] * x_ + p[1]) * x_ + p[2]) * x_ + p[3]);
    return true;
  }

 private:
  const double x_;
  const double y_;
};

// Estimator class.
class CubicEstimator: public theia::Estimator<Point, Cubic> {

public:

  CubicEstimator(){};

  CubicEstimator(int ssize_);

  ~CubicEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Cubic>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Cubic& line) const override;

  private:

  int ssize;

};

/* ****************** */
/* Exponential model. */
/* ****************** */

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

  ExponentialEstimator(int ssize_);

  ~ExponentialEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Exponential>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Exponential& line) const override;

  private:

  int ssize;

};

/* *************** */
/* Logistic model. */
/* *************** */

// Our "model".
struct Logistic {
  double a; double b; double c; double d;
};

// Structure to model LS solver for Ceres. This structure is to find
// linear one-dimensinal models.
struct LogisticResidual {
  LogisticResidual(double x, double y)
      : x_(x), y_(y) {}

  template <typename T> bool operator()(const T* const p,
                                        T* residual) const {
    residual[0] = y_ - (p[0] + p[1] / (1.0 + exp(- p[2] * x_ + p[3])));
    return true;
  }

 private:
  const double x_;
  const double y_;
};

// Estimator class.
class LogisticEstimator: public theia::Estimator<Point, Logistic> {

public:

  LogisticEstimator(){};

  LogisticEstimator(int ssize_);

  ~LogisticEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Logistic>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Logistic& line) const override;

  private:

  int ssize;

};

//////////////////
// Circle model //
//////////////////

// Our "model".
struct Circle {
  double x; double y; double r;
};

// Structure to model LS solver for Ceres. This structure is to find
// linear one-dimensinal models.
struct CircleResidual {
  CircleResidual(double x, double y, double z)
    : x_(x), y_(y), z_(z) {}

  template <typename T> bool operator()(const T* const p,
                                        T* residual) const {
    residual[0] = z_ - (pow(x_ - p[0], 2) + pow(y_ - p[1], 2) - p[2] * p[2]);
    return true;
  }

 private:
  const double x_;
  const double y_;
  const double z_;
};

// Estimator class.
class CircleEstimator: public theia::Estimator<Point, Circle> {

public:

  CircleEstimator() {};
  
  CircleEstimator(int ssize_);

  ~CircleEstimator(){};
  
  // Number of points needed to estimate a line.
  double SampleSize() const override;

  // Estimate a line from two points.
  bool EstimateModel(const std::vector<Point>& data,
                     std::vector<Circle>* models) const override;

  // Calculate the error as the y distance of the point to the line.
  double Error(const Point& point, const Circle& line) const override;

  private:

  int ssize;

};


void run_ransac(void);
