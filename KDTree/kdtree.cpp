// A C++ program to demonstrate operations of KD tree
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
// A structure to represent node of kd tree
class Node
{
public:
    KDPoint point; // To store k dimensional point
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

// A method to create a node of KD tree
Node newNode(KDPoint &point_)
{
    class Node temp;

    temp.point = point_;

    temp.left = temp.right = NULL;
    return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node insertRec(Node *root, KDPoint &point, unsigned depth)
{
    // Tree is empty?
    if (root == nullptr)
        return &newNode(point);

    // Calculate current dimension (cd) of comparison
    unsigned cd = depth % k;

    // Compare the new point with root on current dimension 'cd'
    // and decide the left or right subtree
    if (point.dims[cd] < (root->point.dims[cd]))
        *(root->left) = insertRec(root->left, point, depth + 1);
    else
        *(root->right) = insertRec(root->right, point, depth + 1);

    return *root;
}

// Function to insert a new point with given point in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node insert(Node &root, KDPoint &point)
{
    return insertRec(&root, point, 0);
}

// A utility method to determine if two Points are same
// in K Dimensional space
bool arePointsSame(KDPoint &point1, KDPoint &point2)
{
    // Compare individual pointinate values
    for (int i = 0; i < k; ++i)
        if (point1.dims[i] != point2.dims[i])
            return false;

    return true;
}

// Searches a Point represented by "point[]" in the K D tree.
// The parameter depth is used to determine current axis.
bool searchRec(Node *root, KDPoint &point, unsigned depth)
{
    // Base cases
    if (root == NULL)
        return false;
    if (arePointsSame(root->point, point))
        return true;

    // Current dimension is computed using current depth and total
    // dimensions (k)
    unsigned cd = depth % k;

    // Compare point with root with respect to cd (Current dimension)
    if (point.dims[cd] < root->point.dims[cd])
        return searchRec(root->left, point, depth + 1);

    return searchRec(root->right, point, depth + 1);
}

// Searches a Point in the K D tree. It mainly uses
// searchRec()
bool search(Node *root, KDPoint point)
{
    // Pass current depth as 0
    return searchRec(root, point, 0);
}

// Driver program to test above functions
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
