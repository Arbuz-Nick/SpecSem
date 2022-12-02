// A C++ program to demonstrate operations of KD tree
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
using namespace std;

const int k = 3;

// A structure to represent node of kd tree
struct Node
{
	float point[k]; // To store k dimensional point
	Node *left, *right;
};

// A method to create a node of K D tree
struct Node* newNode(float arr[])
{
	struct Node* temp = new Node;

	for (int i=0; i<k; i++)
	temp->point[i] = arr[i];

	temp->left = temp->right = NULL;
	return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec(Node *root, float point[], unsigned depth)
{
	// Tree is empty?
	if (root == NULL)
	    return newNode(point);

	// Calculate current dimension (cd) of comparison
	unsigned cd = depth % k;

	// Compare the new point with root on current dimension 'cd'
	// and decide the left or right subtree
	if (point[cd] < (root->point[cd]))
		root->left = insertRec(root->left, point, depth + 1);
	else
		root->right = insertRec(root->right, point, depth + 1);

	return root;
}

// Function to insert a new point with given point in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node* insert(Node *root, float point[])
{
	return insertRec(root, point, 0);
}

// A utility method to determine if two Points are same
// in K Dimensional space
bool arePointsSame(float point1[], float point2[])
{
	// Compare individual pointinate values
	for (int i = 0; i < k; ++i)
		if (point1[i] != point2[i])
			return false;

	return true;
}

// Searches a Point represented by "point[]" in the K D tree.
// The parameter depth is used to determine current axis.
bool searchRec(Node* root, float point[], unsigned depth)
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
	if (point[cd] < root->point[cd])
		return searchRec(root->left, point, depth + 1);

	return searchRec(root->right, point, depth + 1);
}

// Searches a Point in the K D tree. It mainly uses
// searchRec()
bool search(Node* root, float point[])
{
	// Pass current depth as 0
	return searchRec(root, point, 0);
}

vector<float> inRectangle(Node* root, float rightTopFront, float leftBottomBack){}

// Driver program to test above functions
int main()
{
	struct Node *root = NULL;
	float points[][k] = {{3, 6, 0}, {17, 15, 0}, {13, 15, 0}, {6, 12, 0},
					{9, 1, 0}, {2, 7, 0}, {10, 19, 0}};

	int n = sizeof(points)/sizeof(points[0]);

	for (int i=0; i<n; i++)
	root = insert(root, points[i]);

	float point1[] = {10, 19, 0};
	(search(root, point1))? cout << "Found\n": cout << "Not Found\n";

	float point2[] = {12, 19, 0};
	(search(root, point2))? cout << "Found\n": cout << "Not Found\n";

	return 0;
}
