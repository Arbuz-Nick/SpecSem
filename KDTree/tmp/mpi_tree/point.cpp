#include "mpi.h"
#include <iostream>
#include <array>
#include <vector>
// #include <opencv2/opencv.hpp>
// #include <opencv2/viz/viz3d.hpp>

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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create the MPI data type for MyClass
    MPI_Datatype myclass_type = MyPoint3D::create_mpi_datatype();

    // Use the new data type to send/receive a MyClass object
    std::vector<MyPoint3D> points(3);
    if (rank == 0)
    {
        MyPoint3D obj(1, 2, 3);
        for (int i = 0; i < 3; i++)
        {
            points[i] = MyPoint3D(i, 1, 2);
        }
        MPI_Send(points.data(), points.size(),
                 myclass_type, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank == 1)
    {
        MyPoint3D recv_obj;
        std::vector<MyPoint3D> rcv_points(3);
        MPI_Recv(rcv_points.data(), rcv_points.size(),
                 myclass_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Do something with the received MyClass object
        for (int i = 0; i < 3; i++)
        {
            std::cout << rcv_points[i][0] << " "
                      << rcv_points[i][1] << " "
                      << rcv_points[i][2] << std::endl;
        }
    }

    MPI_Type_free(&myclass_type);
    MPI_Finalize();
    return 0;
}