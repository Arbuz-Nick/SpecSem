#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <array>
#include <algorithm>
#include <vector>
#include <time.h>

#include "kdtree.h"

#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include <sys/types.h>
#include <signal.h>

#include <opencv2/opencv.hpp>
#include <opencv2/viz/viz3d.hpp>

int axis = 0;
int process_num, rank;
int width = 1000;
int height = 1000;
int deep = 1000;

int side_x = 500;
int side_y = 500;
int side_z = 500;

int p_width = 2;
int l_width = 3;
double bench_t_start, bench_t_end;

static double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday(&Tp, NULL);
    if (stat != 0)
        printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
    bench_t_start = rtclock();
}

void bench_timer_stop()
{
    bench_t_end = rtclock();
}

double bench_timer_print(int rank)
{
    // FILE *fout;
    //  fout = fopen(file_path, "a+");
    //  fprintf(fout, "Rank: %d: %0.6lf;%d;%d*%d*%d*%d\n", rank, bench_t_end - bench_t_start);

    printf("Rank: %d: Time in seconds = %0.10lf\n", rank, bench_t_end - bench_t_start);
    return bench_t_end - bench_t_start;
}

int get_parent(int rank)
{
    return rank / 2 - ((rank + 1) % 2);
}

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

    static MPI_Datatype create_mpi_datatype()
    {
        int block_lengths[] = {3}; // 1 int for a_, 1 double for b_, 10 chars for c_
        MPI_Datatype block_types[] = {MPI_DOUBLE};
        MPI_Aint block_offsets[] = {0};
        MPI_Datatype datatype;
        MPI_Type_create_struct(1, block_lengths, block_offsets, block_types, &datatype);
        MPI_Type_commit(&datatype);
        return datatype;
    }
};
int handler(MyPoint3D lhs, MyPoint3D rhs)
{
    return lhs[axis] < rhs[axis];
}

void MPI_kdtree(MPI_Datatype point_type, std::vector<MyPoint3D> &points, int cloud_size)
{
    int lvl = std::log2(rank + 1);

    int npoints = 0;
    int parent_rank = get_parent(rank);
    int left_child = rank * 2 + 1;
    int right_child = rank * 2 + 2;
    if (left_child >= process_num)
    {
        left_child = -1;
        right_child = -1;
    }
    else if (right_child >= process_num)
    {
        right_child = -1;
    }

    axis = lvl % MyPoint3D::DIM;

    if (rank == 0)
    {
        for (int i = 0; i < cloud_size; i++)
        {
            const int x = rand() % width;
            const int y = rand() % height;
            const int z = rand() % deep;
            points.push_back(MyPoint3D(x, y, z));
        }
        npoints = points.size();

        std::cout << "Main is sending\n";
        if (right_child != -1)
        {

            const int mid = (npoints - 1) / 2;
            std::nth_element(points.data(), points.data() + mid, points.data() + npoints,
                             handler);
            const int left_size = npoints / 2;
            ;
            const int right_size = npoints - left_size;

            MPI_Send(&left_size, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data(), left_size, point_type, left_child, 1, MPI_COMM_WORLD);

            MPI_Send(&right_size, 1, MPI_INT, right_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data() + mid, right_size, point_type, right_child, 1, MPI_COMM_WORLD);
        }
        else if (left_child != -1)
        {
            MPI_Send(&npoints, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data(), npoints, point_type, left_child, 1, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&npoints, 1, MPI_INT, parent_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < npoints; i++)
        {
            points.push_back(MyPoint3D(0, 0, 0));
        }
        MPI_Recv(points.data(), npoints, point_type, parent_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        npoints = points.size();

        if (right_child != -1)
        {

            const int mid = (npoints - 1) / 2;

            std::nth_element(points.data(), points.data() + mid, points.data() + npoints,
                             handler);

            const int left_size = npoints / 2;
            ;
            const int right_size = npoints - left_size;

            MPI_Send(&left_size, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data(), left_size, point_type, left_child, 1, MPI_COMM_WORLD);

            MPI_Send(&right_size, 1, MPI_INT, right_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data() + mid, right_size, point_type, right_child, 1, MPI_COMM_WORLD);
        }
        else if (left_child != -1)
        {
            MPI_Send(&npoints, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD);
            MPI_Send(points.data(), npoints, point_type, left_child, 1, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype point_type = MyPoint3D::create_mpi_datatype();
    int lvl = std::log2(rank + 1);

    int parent_rank = get_parent(rank);
    int left_child = rank * 2 + 1;
    int right_child = rank * 2 + 2;
    if (left_child >= process_num)
    {
        left_child = -1;
        right_child = -1;
    }
    else if (right_child >= process_num)
    {
        right_child = -1;
    }

    std::vector<MyPoint3D> points;
    if (rank == 0)
        bench_timer_start();
    MPI_kdtree(point_type, points, 10000);
    if (rank == 0)
    {
        bench_timer_stop();
        bench_timer_print(rank);
    }
    if (lvl == int(std::log2(process_num)))
    {
        bench_timer_start();
        std::cout << "My rank: " << rank << ", Got from " << parent_rank << ", data size: " << points.size() << std::endl;
        kdt::KDTree<MyPoint3D> kdtree(points);
        /*cv::viz::Viz3d window;
        std::vector<cv::Point3d> points_3d_(points.size());
        for (int i = 0; i < points.size(); i++)
            points_3d_[i] = points[i];
        if (points_3d_.empty())
            std::cout << "My rank: " << rank << "Cloud is empty\n";
        cv::viz::WCloud point_cloud = cv::viz::WCloud(points_3d_, cv::viz::Color::green());
        point_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
        window.showWidget("points", point_cloud);
        */
        const MyPoint3D bot_3d((width - side_x) / 2, (height - side_y) / 2, (deep - side_z) / 2),
            top_3d((width + side_x) / 2, (height + side_y) / 2, (deep + side_z) / 2);
        /*
        cv::viz::WCube cube = cv::viz::WCube(bot_3d, top_3d, true, cv::viz::Color::blue());
        cube.setRenderingProperty(cv::viz::LINE_WIDTH, l_width);
        window.showWidget("cube", cube);
        */
        const std::vector<int> cubeIndices = kdtree.rectangleSearch(bot_3d, top_3d);
        std::vector<MyPoint3D> found_points_3d;
        for (auto i : cubeIndices)
        {
            found_points_3d.push_back(points[i]);
        }
        /*
        if (!found_points_3d.empty())
        {
            cv::viz::WCloud found_cloud = cv::viz::WCloud(found_points_3d, cv::viz::Color::red());
            found_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
            window.showWidget("found_points", found_cloud);
        }
        */
        // std::cout << points_3d[10][0] << std::endl;

        bench_timer_stop();
        // std::cout << "My rank: " << rank << ", Iter per point: " << double(kdt::niter) / points.size() << std::endl;
        double t = bench_t_end - bench_t_start;
        MPI_Send(&t, 1, MPI_DOUBLE, parent_rank, 1, MPI_COMM_WORLD);
        int i_size = found_points_3d.size();
        MPI_Send(&i_size, 1, MPI_INT, parent_rank, 1, MPI_COMM_WORLD);
        MPI_Send(found_points_3d.data(), i_size, point_type, parent_rank, 1, MPI_COMM_WORLD);

        // window.spin();
    }
    else
    {
        // Получение времени
        double r_t = 0, l_t = 0;
        if (right_child != -1)
        {
            MPI_Recv(&l_t, 1, MPI_DOUBLE, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&r_t, 1, MPI_DOUBLE, right_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Recv(&l_t, 1, MPI_DOUBLE, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        double max_t = std::max(l_t, r_t);
        if (parent_rank != -1)
            MPI_Send(&max_t, 1, MPI_DOUBLE, parent_rank, 1, MPI_COMM_WORLD);
        else
            printf("Rank: %d: Max time in seconds = %0.10lf\n", rank, max_t);

        // Получение точек
        points.clear();
        points.shrink_to_fit();

        int r_size = 0, l_size = 0;
        if (right_child != -1)
        {
            MPI_Recv(&l_size, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < l_size; i++)
                points.push_back(MyPoint3D());
            MPI_Recv(points.data(), l_size, point_type, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Recv(&r_size, 1, MPI_INT, right_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < r_size; i++)
                points.push_back(MyPoint3D());
            MPI_Recv(points.data() + l_size, r_size, point_type, right_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Recv(&l_size, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < l_size; i++)
                points.push_back(MyPoint3D());
            MPI_Recv(points.data(), l_size, point_type, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (parent_rank != -1)
        {
            int size = points.size();
            MPI_Send(&size, 1, MPI_INT, parent_rank, 1, MPI_COMM_WORLD);
            MPI_Send(points.data(), size, point_type, parent_rank, 1, MPI_COMM_WORLD);
        }
        else
        {
            // printf("Rank: %d: Max time in seconds = %0.10lf\n", rank, max_t);
            cv::viz::Viz3d window;
            std::vector<cv::Point3d> points_3d_(points.size());
            for (int i = 0; i < points.size(); i++)
                points_3d_[i] = points[i];
            if (points_3d_.empty())
                std::cout << "My rank: " << rank << "Cloud is empty\n";
            cv::viz::WCloud point_cloud = cv::viz::WCloud(points_3d_, cv::viz::Color::green());
            point_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
            window.showWidget("points", point_cloud);

            const MyPoint3D bot_3d((width - side_x) / 2, (height - side_y) / 2, (deep - side_z) / 2), top_3d((width + side_x) / 2, (height + side_y) / 2, (deep + side_z) / 2);
            cv::viz::WCube cube = cv::viz::WCube(bot_3d, top_3d, true, cv::viz::Color::blue());
            cube.setRenderingProperty(cv::viz::LINE_WIDTH, l_width);
            window.showWidget("cube", cube);

            //std::vector<cv::Point3d> found_points_3d;
            /*for (auto i : indexes)
            {
                found_points_3d.push_back(points_3d_[i]);
            }
            if (!found_points_3d.empty())
            {
                cv::viz::WCloud found_cloud = cv::viz::WCloud(found_points_3d, cv::viz::Color::red());
                found_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
                window.showWidget("found_points", found_cloud);
            }*/
            window.spin();
        }
    }
    /*
    if (rank == 0)
    {
        bench_timer_start();
        std::cout << "My rank: " << rank << ", Got from " << parent_rank << ", data size: " << points.size() << std::endl;
        kdt::KDTree<MyPoint3D> kdtree(points);
        cv::viz::Viz3d window;
        std::vector<cv::Point3d> points_3d_(points.size());
        for (int i = 0; i < points.size(); i++)
            points_3d_[i] = points[i];
        if (points_3d_.empty())
            std::cout << "My rank: " << rank << "Cloud is empty\n";
        cv::viz::WCloud point_cloud = cv::viz::WCloud(points_3d_, cv::viz::Color::green());
        point_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
        window.showWidget("points", point_cloud);

        const MyPoint3D bot_3d((width - side_x) / 2, (height - side_y) / 2, (deep - side_z) / 2), top_3d((width + side_x) / 2, (height + side_y) / 2, (deep + side_z) / 2);
        cv::viz::WCube cube = cv::viz::WCube(bot_3d, top_3d, true, cv::viz::Color::blue());
        cube.setRenderingProperty(cv::viz::LINE_WIDTH, l_width);
        window.showWidget("cube", cube);

        const std::vector<int> cubeIndices = kdtree.rectangleSearch(bot_3d, top_3d);
        std::vector<cv::Point3d> found_points_3d;
        for (auto i : cubeIndices)
        {
            found_points_3d.push_back(points_3d_[i]);
        }
        if (!found_points_3d.empty())
        {
            cv::viz::WCloud found_cloud = cv::viz::WCloud(found_points_3d, cv::viz::Color::red());
            found_cloud.setRenderingProperty(cv::viz::POINT_SIZE, p_width);
            window.showWidget("found_points", found_cloud);
        }
        // std::cout << points_3d[10][0] << std::endl;
        bench_timer_stop();
        std::cout << "My rank: " << rank << ", Iter per point: " << double(kdt::niter) / points.size() << std::endl;
        double t = bench_timer_print(rank);
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(&point_type);
    MPI_Finalize();
    return 0;
}

/*
               0
        1              2
   3       4       5       6
 7   8   9  10  11  12  13  14
*/