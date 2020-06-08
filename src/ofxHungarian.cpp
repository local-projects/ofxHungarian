#include "ofxHungarian.h"

// --------------------------------------------------
ofxHungarian::ofxHungarian() {

}

// --------------------------------------------------
ofxHungarian::~ofxHungarian() {
}

// --------------------------------------------------
void ofxHungarian::solve(vector<HungarianSample>& from, vector<HungarianSample>& to, float radius) {

	// The type of our cost function is float
	typedef float Cost;

	// Create an array to store the costs
	Matrix2D<Cost> cost_matrix(to.size(), from.size(), 0);

	// Calculate the amount of intersection over each point with a bounding box with given radius 
	// on each side.
	for (int row = 0; row < to.size(); row++) {
		for (int col = 0; col < from.size(); col++) {
			// invert the intersection so we minimize cost
			cost_matrix(row, col) = -IOU(to[row].position, from[col].position, radius);
		}
	}

	// Create an assignment matrix for all of the samples we're mapping to
	vector<size_t> assignment(cost_matrix.rows());

	// Solve
	//uint64_t startTime = ofGetElapsedTimeMicros();
	Cost total_cost = hungarian_algorithm::solve<Cost>(cost_matrix, cost_matrix.rows(),
		cost_matrix.cols(), assignment);
	//uint64_t elapsedTime = ofGetElapsedTimeMicros() - startTime;
	//cout << "Took " << float(elapsedTime) / 1000.0 << " ms" << endl;

	// Iterate through all assignments.
	// If valid, save the link
	for (int row = 0; row < assignment.size(); row++) {

		int col = assignment[row];
		if (col >= 0 && col < from.size() && cost_matrix(row, col) < 0) {
			// Valid assignment; Save the index in the from array.
			to[row].mapTo = col;
		}
	}
}

// --------------------------------------------------
float ofxHungarian::IOU(glm::vec3& a, glm::vec3& b, float r) {

	// Check if there exists an intersection. If not, return false.
	for (int i = 0; i < 3; i++) {
		if ((a[i] - r) > (b[i] + r) || (a[i] + r) < (b[i] - r)) return 0;
	}
	// Alt: if distance between points is greater than r
	float iou = (2 * r - abs(a.x - b.x))* (2 * r - abs(a.y - b.y))* (2 * r - abs(a.z - b.z)) / pow(2 * r, 3);
	return iou;
}

// --------------------------------------------------

// --------------------------------------------------

// --------------------------------------------------

// --------------------------------------------------

