#include <iostream>
#include <cmath>
#include <vector>
#include <sys/types.h>
#include <unistd.h>
#include <numeric>
#include <algorithm>
#include <functional>
int handler(int lhs, int rhs)
{
    return lhs < rhs;
}
int main(int argc, char const *argv[])
{
    int level, rank, n = 63;
    // for (int i = 0; i < n; i++)
    int n_level = std::log2(n);

    // printf("%d\n", int(std::log2(1020)));
    std::vector<int> points(n);
    for (int i = 0; i < n; i++)
    {
        points[i] = 30 - i + 1;
    }
    const int mid = (n - 1) / 2;
    std::vector<int> indices(points.size());
    std::iota(std::begin(indices), std::end(indices), 0);

    std::nth_element(
        points.begin(), points.begin() + mid, points.end(),
        handler);
    {
        //   return lhs < rhs;
    };
    for (int i = 0; i < n; i++)
    {
        std::cout << points[i] << std::endl;
    }
    return 0;
}
