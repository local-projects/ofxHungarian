#pragma once

#include "ofMain.h"
#include "hungarian_algorithm.h"
#include "Matrix2D.h"

class HungarianSample {
public:

    string key = "";
    glm::vec3 position;
    int index;

    // Will this sample be mapped to another sample?
    // If so, this index should be >= 0
    int mapTo = -1;

};

// A filter manipulates realtime data using a series of ops (operations)
class ofxHungarian {
public:
    
    ofxHungarian();
    ~ofxHungarian();

    // This solver is designed to solve for point samples, using the radius provided
    // as a common measure of proximity.
    static void solve(vector<HungarianSample>& from, vector<HungarianSample>& to, float radius);

private:

    // Calculate intersection over union.
    // Two points a and b with boxes with radii r.
    static float IOU(glm::vec3& a, glm::vec3& b, float r);

};
