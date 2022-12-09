#include <opencv2/opencv.hpp>
#include <opencv2/viz/viz3d.hpp>
//#include <opencv2/viz/viz3d.hpp>
//#include <opencv2/viz/vizcore.hpp>

#include <iostream>
#include <array>
#include <vector>

#include "kdtree.h"

class MyPoint : public std::array<double, 2>
{
public:

	// dimension of space (or "k" of k-d tree)
	// KDTree class accesses this member
	static const int DIM = 2;

	// the constructors
	MyPoint() {}
	MyPoint(double x, double y)
	{ 
		(*this)[0] = x;
		(*this)[1] = y;
	}

	// conversion to OpenCV Point2d
	operator cv::Point2d() const { return cv::Point2d((*this)[0], (*this)[1]); }
};

// user-defined point type
// inherits std::array in order to use operator[]
class MyPoint3D : public std::array<double, 3>
{
public:

	// dimension of space (or "k" of k-d tree)
	// KDTree class accesses this member
	static const int DIM = 3;
    
	// the constructors
	MyPoint3D() {}
	MyPoint3D(double x, double y, double z)
	{ 
		(*this)[0] = x;
		(*this)[1] = y;
        (*this)[2] = z;
	}

	// conversion to OpenCV Point2d
	operator cv::Point3d() const { return cv::Point3d((*this)[0], (*this)[1], (*this)[2]); }
};

int main(int argc, char **argv)
{

	const int seed = argc > 1 ? std::stoi(argv[1]) : 0;
	srand(seed);

	// generate space
	int width = 500;
	int height = 500;
	cv::Mat img = cv::Mat::zeros(cv::Size(width, height), CV_8UC3);

	// generate points
	int npoints = 100;
	std::vector<MyPoint> points(npoints);
	for (int i = 0; i < npoints; i++)
	{
		const int x = rand() % width;
		const int y = rand() % height;
		points[i] = MyPoint(x, y);
	}

	for (const auto& pt : points)
		cv::circle(img, cv::Point2d(pt), 1, cv::Scalar(0, 255, 255), -1);

	// build k-d tree
	kdt::KDTree<MyPoint> kdtree(points);

	// generate query (center of the space)
	const MyPoint query(0.5 * width, 0.5 * height);
	cv::circle(img, cv::Point2d(query), 1, cv::Scalar(0, 0, 255), -1);

	// nearest neigbor search
	const cv::Mat I0 = img.clone();
	const int idx = kdtree.nnSearch(query);
	cv::circle(I0, cv::Point2d(points[idx]), 1, cv::Scalar(255, 255, 0), -1);
	cv::line(I0, cv::Point2d(query), cv::Point2d(points[idx]), cv::Scalar(0, 0, 255));

	// k-nearest neigbors search
	const cv::Mat I1 = img.clone();
	const int k = 10;
	const std::vector<int> knnIndices = kdtree.knnSearch(query, k);
	for (int i : knnIndices)
	{
		cv::circle(I1, cv::Point2d(points[i]), 1, cv::Scalar(255, 255, 0), -1);
		cv::line(I1, cv::Point2d(query), cv::Point2d(points[i]), cv::Scalar(0, 0, 255));
	}
	
	// radius search
	const cv::Mat I2 = img.clone();
	const double radius = 50;
	const std::vector<int> radIndices = kdtree.radiusSearch(query, radius);
	for (int i : radIndices)
		cv::circle(I2, cv::Point2d(points[i]), 1, cv::Scalar(255, 255, 0), -1);
	cv::circle(I2, cv::Point2d(query), cvRound(radius), cv::Scalar(0, 0, 255));

    const cv::Mat I3 = img.clone();
    int side = 200;
    const MyPoint bot(width*0.5 - side/2, height*0.5 - side/2), top(width*0.5 + side/2, height*0.5 + side/2);
    const std::vector<int> rectIndices = kdtree.rectangleSearch(bot, top);
    cv::Rect rect(width*0.5 - side/2, height*0.5 - side/2, side, side);
	for (int i : rectIndices)
		cv::circle(I3, cv::Point2d(points[i]), 1, cv::Scalar(255, 255, 0), -1);

    
	cv::rectangle(I3, rect, cv::Scalar(0, 0, 255));


	int side_x = 100;
    int side_y = 200;
    int side_z = 100;

    width = 1000;
    height = 1000;
    int deep = 1000;
    npoints *= 10;
    int p_width = 2;
    int l_width = 3;
    

    //cv::viz::Viz3d window; //creating a Viz window
	//Displaying the Coordinate Origin (0,0,0)
	//window.showWidget("coordinate", cv::viz::WCoordinateSystem(100));
	//Displaying the 3D points in green
    for (int npoints_i = 10; npoints_i <= 100000; npoints_i *= 10){
        kdt::niter = 0;
        cv::viz::Viz3d window;
        std::vector<MyPoint3D> points_3d(npoints_i);
        for (int i = 0; i < npoints_i; i++)
        {
            const int x = rand() % width;
            const int y = rand() % height;
            const int z = rand() % deep;
            points_3d[i] = MyPoint3D(x, y, z);
        }
        std::vector<cv::Point3d> points_3d_(npoints_i);
        for (int i = 0; i < npoints_i; i++)
            points_3d_[i] = points_3d[i];
        kdt::KDTree<MyPoint3D> kdtree_3d(points_3d);
        cv::viz::WCloud point_cloud = cv::viz::WCloud(points_3d_, cv::viz::Color::green());
        point_cloud.setRenderingProperty( cv::viz::POINT_SIZE, p_width );
        window.showWidget("points", point_cloud);

        const MyPoint3D bot_3d((width - side_x)/2, (height - side_y)/2, (deep - side_z)/2), top_3d((width + side_x)/2, (height + side_y)/2, (deep + side_z)/2);
        cv::viz::WCube cube = cv::viz::WCube(bot_3d, top_3d, true, cv::viz::Color::blue());
        cube.setRenderingProperty( cv::viz::LINE_WIDTH, l_width );
        window.showWidget("cube", cube);
        
        const std::vector<int> cubeIndices = kdtree_3d.rectangleSearch(bot_3d, top_3d);
        std::vector<cv::Point3d> found_points_3d;
        for (auto i : cubeIndices) {
            found_points_3d.push_back(points_3d_[i]);
        }
        cv::viz::WCloud found_cloud = cv::viz::WCloud(found_points_3d, cv::viz::Color::red());
        found_cloud.setRenderingProperty( cv::viz::POINT_SIZE, p_width );
        window.showWidget("found_points", found_cloud);
        //std::cout << points_3d[10][0] << std::endl;
        std::cout << double(kdt::niter) / npoints_i << std::endl;
        window.spin();
    }

	// show results

	cv::imshow("Nearest neigbor search", I0);
	cv::imshow("K-nearest neigbors search (k = 10)", I1);
	cv::imshow("Radius search (radius = 50)", I2);
	cv::imshow("Rectangle search (side = 100)", I3);


	cv::waitKey(0);

	return 0;
}