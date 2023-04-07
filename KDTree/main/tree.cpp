#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <array>
#include <algorithm>
#include <vector>
#include <time.h>
#include <string>

#include "kdtree.h"

#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include <sys/types.h>
#include <signal.h>


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
double left_border[2], right_border[2];
double build_time, search_time;

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

    static MPI_Datatype create_mpi_datatype()
    {
        int block_lengths[] = {3};
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

void MPI_kdtree(MPI_Datatype point_type, std::vector<MyPoint3D> &points, int cloud_size = -1, std::string path = "")
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
        if (path == "")
            if (cloud_size == -1)
            {
                std::cout << "ERROR:: Can't build cloud of such size" << std::endl;
                return;
            }
            else
                for (int i = 0; i < cloud_size; i++)
                {
                    const int x = rand() % width;
                    const int y = rand() % height;
                    const int z = rand() % deep;
                    points.push_back(MyPoint3D(x, y, z));
                }
        else
        {
            double x, y, z;
            std::ifstream file;
            file.open(path);
            if (!file.is_open())
            {
                std::cout << "File opening error\n"
                          << path << std::endl;
                return;
            }
            while (!file.eof())
            {
                file >> x >> y >> z;
                points.push_back(MyPoint3D(x, y, z));
            }
        }

        npoints = points.size();

        if (right_child != -1)
        {

            const int mid = (npoints - 1) / 2;
            std::nth_element(points.data(), points.data() + mid, points.data() + npoints,
                             handler);
            std::nth_element(points.data(), points.data(), points.data() + npoints,
                             handler);
            std::nth_element(points.data(), points.data() + npoints - 1, points.data() + npoints,
                             handler);

            const int left_size = npoints / 2;
            const int right_size = npoints - left_size;

            left_border[0] = points[0][axis];
            left_border[1] = points[mid][axis];
            right_border[0] = points[mid][axis];
            right_border[1] = points[npoints - 1][axis];

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
            std::nth_element(points.data(), points.data(), points.data() + npoints,
                             handler);
            std::nth_element(points.data(), points.data() + npoints - 1, points.data() + npoints,
                             handler);

            const int left_size = npoints / 2;

            const int right_size = npoints - left_size;
            left_border[0] = points[0][axis];
            left_border[1] = points[mid][axis];
            right_border[0] = points[mid][axis];
            right_border[1] = points[npoints - 1][axis];

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
    if ((left_child != -1) && (rank != 0))
    {
        points.clear();
        points.shrink_to_fit();
    }
}

void rectangleSearch(MPI_Datatype point_type, MyPoint3D borders[2], std::vector<MyPoint3D> &points)
{
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
    
    if (left_child == -1)
    {

        if (parent_rank != -1)
        {
            MPI_Recv(borders, 2, point_type, parent_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if ((borders[0][0] == borders[1][0]) && (borders[0][1] == borders[1][1]) && (borders[0][2] == borders[1][2]))
        {
            return;
        }
        kdt::KDTree<MyPoint3D> kdtree(points);
        
        const std::vector<int> cubeIndices = kdtree.rectangleSearch(borders[0], borders[1]);
        std::vector<MyPoint3D> found_points_3d;
        for (auto i : cubeIndices)
        {
            found_points_3d.push_back(points[i]);
        }
        if (parent_rank != -1)
        {
            int i_size = found_points_3d.size();
            MPI_Send(&i_size, 1, MPI_INT, parent_rank, 1, MPI_COMM_WORLD);
            MPI_Send(found_points_3d.data(), i_size, point_type, parent_rank, 1, MPI_COMM_WORLD);
        }
        else
        {
            points = found_points_3d;
        }
    }
    else
    {

        if (rank != 0)
        {
            MPI_Recv(borders, 2, point_type, parent_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (right_child != -1)
        {
            if (borders[0][axis] <= left_border[1])
            {
                MPI_Send(borders, 2, point_type, left_child, 1, MPI_COMM_WORLD);
            }
            else
            {
                MyPoint3D not_borders[2];
                not_borders[0] = borders[1];
                not_borders[1] = borders[1];
                MPI_Send(not_borders, 2, point_type, left_child, 1, MPI_COMM_WORLD);
            }
            if (borders[1][axis] >= right_border[0])
            {
                MPI_Send(borders, 2, point_type, right_child, 1, MPI_COMM_WORLD);
            }
            else
            {
                MyPoint3D not_borders[2];
                not_borders[0] = borders[1];
                not_borders[1] = borders[1];
                MPI_Send(not_borders, 2, point_type, right_child, 1, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Send(borders, 2, point_type, left_child, 1, MPI_COMM_WORLD);
        }

        if ((borders[0][0] == borders[1][0]) && (borders[0][1] == borders[1][1]) && (borders[0][2] == borders[1][2]))
        {
            return;
        }

        // Получение точек
        int r_size = 0, l_size = 0;
        if (right_child != -1)
        {
            if (borders[0][axis] <= left_border[1])
            {
                MPI_Recv(&l_size, 1, MPI_INT, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < l_size; i++)
                    points.push_back(MyPoint3D());
                MPI_Recv(points.data(), l_size, point_type, left_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (borders[1][axis] >= right_border[0])
            {
                MPI_Recv(&r_size, 1, MPI_INT, right_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < r_size; i++)
                    points.push_back(MyPoint3D());
                MPI_Recv(points.data() + l_size, r_size, point_type, right_child, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
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
    }
}

int main(int argc, char *argv[])
{
    int seed;
    sscanf(argv[2], "%d", &seed);
    std::srand(time(NULL));

    int point_number = -1;
    std::string path = "";

    switch (argv[1][0])
    {
    case 'd':
        path = "./../data/Fantasy Dragon.xyz";
        break;
    case 'f':
        path = "./../data/Bearded guy.xyz";
        break;
    default:
        sscanf(argv[1], "%d", &point_number);
        break;
    }

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype point_type = MyPoint3D::create_mpi_datatype();

    std::vector<MyPoint3D> points;
    if (rank == 0)
        bench_timer_start();

    MPI_kdtree(point_type, points, point_number, path);
    MPI_Barrier(MPI_COMM_WORLD);
    MyPoint3D borders[2];
    if (rank == 0)
    {
        bench_timer_stop();
        build_time = bench_t_end - bench_t_start;
        bench_timer_start();
        MyPoint3D p = points[std::rand() % points.size()];
        points.clear();
        points.shrink_to_fit();
        borders[0] = p;
        borders[1] = MyPoint3D(p[0] + std::rand() % 500, p[1] + std::rand() % 500, p[2] + std::rand() % 500);
    }

    rectangleSearch(point_type, borders, points);

    if (rank == 0)
    {
        bench_timer_stop();
        search_time = bench_t_end - bench_t_start;
        std::ofstream out_file;
        out_file.open("result_mpi_polus.csv", std::ios_base::app);
        out_file << build_time << ";" << search_time << ";" << process_num << ";" << argv[1] << ";" << points.size() << std::endl;
        out_file.close();

    }

    MPI_Type_free(&point_type);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
