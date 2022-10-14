#include <vector>
#include <iostream>

const int k = 3;
class KDPoint
{
public:
    std::vector<float> dims;
    KDPoint(float x_ = 0, float y_ = 0, float z_ = 0)
    {
        dims.push_back(x_);
        dims.push_back(y_);
        dims.push_back(z_);
    }
    KDPoint(std::vector<float> vector_)
    {
        for (int i = 0; i < k; i++)
        {
            dims.push_back(vector_[i]);
        }
    }
    KDPoint(KDPoint &point_)
    {
        dims = point_.dims;
    }

    KDPoint &operator=(KDPoint const &point_)
    {
        if (this == &point_)
            return *this;
        this->dims = point_.dims;
        return *this;
    }
};
class Node
{
public:
    KDPoint point;
    Node *left, *right;
    Node()
    {
        left = nullptr;
        right = nullptr;
    }
    Node(Node *node_)
    {
        left = node_->left;
        right = node_->right;
        point = node_->point;
    }
    Node operator=(Node const &node_)
    {
        if (this == &node_)
            return *this;
        this->point = node_.point;
        this->left = node_.left;
        this->right = node_.right;
        return *this;
    }
};

Node newNode(KDPoint &point_)
{
    class Node temp;

    temp.point = point_;

    temp.left = temp.right = NULL;
    return temp;
}

Node insertRec(Node *root, KDPoint &point, unsigned depth)
{
    if (root == nullptr)
        return &newNode(point);

    unsigned cd = depth % k;

    if (point.dims[cd] < (root->point.dims[cd]))
        *(root->left) = insertRec(root->left, point, depth + 1);
    else
        *(root->right) = insertRec(root->right, point, depth + 1);

    return *root;
}

Node insert(Node &root, KDPoint &point)
{
    return insertRec(&root, point, 0);
}

bool arePointsSame(KDPoint &point1, KDPoint &point2)
{
    for (int i = 0; i < k; ++i)
        if (point1.dims[i] != point2.dims[i])
            return false;

    return true;
}

bool searchRec(Node *root, KDPoint &point, unsigned depth)
{
    if (root == NULL)
        return false;
    if (arePointsSame(root->point, point))
        return true;

    unsigned cd = depth % k;

    if (point.dims[cd] < root->point.dims[cd])
        return searchRec(root->left, point, depth + 1);

    return searchRec(root->right, point, depth + 1);
}

bool search(Node *root, KDPoint point)
{
    return searchRec(root, point, 0);
}

int main()
{
    struct Node *root = NULL;
    int points[][k] = {{3, 6, 1}, {17, 15, 2}, {13, 15, 3}, {6, 12, 4}, {9, 1, 5}, {2, 7, 6}, {10, 19, 7}};
    std::vector<KDPoint> point_vector;
    for (int i = 0; i < sizeof(points) / sizeof(points[0]); ++i)
    {
        point_vector.emplace_back(point_vector[i, 0], point_vector[i, 1], point_vector[i, 2]);
    }
    for (auto point : point_vector)
    {
        root = &insert(*root, point);
    }

    KDPoint point1(10, 19, 7);
    (search(root, point1)) ? std::cout << "Found\n" : std::cout << "Not Found\n";

    KDPoint point2(12, 19, 0);
    (search(root, point2)) ? std::cout << "Found\n" : std::cout << "Not Found\n";

    return 0;
}
