# To build and compile Theia
#
# Create a cmake/FindCeres.cmake with the following lines
#
# SET(CERES_FOUND TRUE)
# SET(CERES_LIBRARY "/usr/local/lib/libceres.so")
# SET(CERES_LIBRARIES "/usr/local/lib/libceres.so")
# SET(CERES_INCLUDE_DIRS "/usr/local/include/ceres/")
#
# Create cmake/FindRapidJSON.cmake
#
# SET(RapidJSON_FOUND TRUE)
# SET(RAPIDJSON_INCLUDE_DIRS "/usr/include/rapidjson")
#
# Run
#
# cmake -DSTATX_WITH_CERES=ON -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS="-shared -fPIC" ..
#
# Edit
#
# applications/CMakeFiles/view_reconstruction.dir/link.txt
#
# and remove the -lFALSE flags. This is the shared library version,
# therefore the Ceres library should also be compiled as a shared
# library, but it is easier to compile.
#
# There are other minor issues that can be easily solved.

using Cxx
using Libdl

push!(DL_LOAD_PATH, "/usr/local/lib")

function ransac_example()

    for header in ["/usr/include/eigen3/",
                   "/usr/local/include/theia/libraries/vlfeat",
                   "/usr/include/suitesparse",
                   "/usr/local/include/theia/libraries/",
                   "/usr/local/include/theia/libraries/statx",
                   "/usr/local/include/theia/libraries/optimo",
                   "/usr/local/include",
                   "/usr/local/include/theia",
                   "/usr/include"]
    
        addHeaderDir(header, kind=C_User)

    end

    Libdl.dlopen("libvlfeat.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libflann_cpp.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libakaze.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libceres.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libtheia.so", Libdl.RTLD_GLOBAL)

    cxx"""
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
    """
    
end

"""

Test the creating and call of a C++ function in Julia. Should print
10.

"""
function hello_world()

    cxx""" #include<iostream> """
    
    cxx"""
           void printme(int x) {
              std::cout << x << std::endl;
           }
    """

    @cxx printme(10)

end
