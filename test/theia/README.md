How to install Theia
====================

We are using Theia v0.8.

  - Create a cmake/FindCeres.cmake with the following lines

        SET(CERES_FOUND TRUE)
        SET(CERES_LIBRARY "/usr/local/lib/libceres.so")
        SET(CERES_LIBRARIES "/usr/local/lib/libceres.so")
        SET(CERES_INCLUDE_DIRS "/usr/local/include/ceres/")

  - Create cmake/FindRapidJSON.cmake

        SET(RapidJSON_FOUND TRUE)
        SET(RAPIDJSON_INCLUDE_DIRS "/usr/include/rapidjson")


  - Install necessary libraries:

        sudo aptitude install libopenimageio-dev libjpeg-dev libgl1-mesa-dev freeglut3-dev libflann-dev
        
  - Apply the fix suggested in [this pull request](https://github.com/sweeneychris/TheiaSfM/pull/227) in source file.
  
  - Run commands
  
        mkdir build
        cmake -DSTATX_WITH_CERES=ON -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS="-shared -fPIC" ..
        make -j4
        
  - Fix linking file `applications/CMakeFiles/view_reconstruction.dir/link.txt` in the `build` directory and remove all `-lFALSE` flags.

  - Install package
  
        sudo make install

How to compile Theia
====================

Run command

    g++ -std=c++11 -I/usr/include/eigen3 -I/usr/local/include/theia/libraries/vlfeat -I/usr/include/suitesparse -I/usr/local/include/theia/libraries/ -I/usr/local/include/theia/libraries/statx -I/usr/local/include/theia/libraries/optimo -I/usr/local/include/theia PROGRAM.CPP -o PROGRAM -ltheia -lglog -lceres
