#include "mpi.h"
#include <iostream>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <array>
#include <algorithm>
#include <vector>
int axis = 0;
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
    // operator cv::Point3d() const { return cv::Point3d((*this)[0], (*this)[1], (*this)[2]); }

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
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int process_num, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype point_type = MyPoint3D::create_mpi_datatype();

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

    int lvl = std::log2(rank + 1);

    int send_buf[lvl + 1];
    send_buf[lvl] = rank;

    // std::cout << "My rank: " << rank << ", Left = " << left_child << ", Right = " << right_child << ", PID = " << getpid() << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    int npoints = 0;
    if (rank == 0)
    {
        std::vector<MyPoint3D> points;
        for (int i = 0; i < 1500; i++)
        {
            points.push_back(MyPoint3D(i, 1, 2));
        }
        axis = lvl % MyPoint3D::DIM;
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
        std::vector<MyPoint3D> points(npoints);
        MPI_Recv(points.data(), npoints, point_type, parent_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (lvl == int(std::log2(process_num)))
        {
            std::cout << "My rank: " << rank << ", Got from " << parent_rank << ", data size: " << points.size() << std::endl;
        }
        axis = lvl % MyPoint3D::DIM;
        npoints = points.size();

        std::cout << "Not Main is sending\n";
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