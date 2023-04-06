#include <fstream>
#include <iostream>
#include <vector>

#include <opencv2/opencv.hpp>
#include <opencv2/viz/viz3d.hpp>
int main(int argc, char const *argv[])
{
    double x, y, z;
    std::ifstream file;
    file.open(argv[1]);
    if (!file.is_open())
    {
        std::cout << "File opening error\n"
                  << argv[1] << std::endl;
        return 1;
    }
    std::cout << "Start to read the file" << std::endl;
    std::vector<cv::Point3d> points_3d;
    // for (int i = 0; i < points.size(); i++)
    //    points_3d[i] = points[i];
    cv::viz::Viz3d window("Show");
    int i = 0;
    while (!file.eof())
    {
        file >> x >> y >> z;
        if (i % 2)
        {
            points_3d.push_back(cv::Point3d(x, y, z));
        }
        i += 1;
    }
    cv::viz::WCloud point_cloud = cv::viz::WCloud(points_3d, cv::viz::Color::green());
    point_cloud.setRenderingProperty(cv::viz::POINT_SIZE, 1);
    window.showWidget("points", point_cloud);
    window.spin();
    return 0;
}
