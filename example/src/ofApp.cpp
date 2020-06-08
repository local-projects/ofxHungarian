#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){

    //int size = 100;

   //typedef float Cost;

   //Matrix2D<Cost> cost_matrix(size, size);
   //for (int i = 0; i < size; i++) {
   //    for (int j = 0; j < size; j++) {
   //        cost_matrix(i, j) = Cost(ofRandom(size * size));
   //    }
   //}
   ////cost_matrix(0, 0) = 1;
   ////cost_matrix(0, 1) = 0;
   ////cost_matrix(1, 0) = 1;
   ////cost_matrix(1, 1) = 2;
   //cout << cost_matrix << endl;

   //vector<size_t> assignment(cost_matrix.rows());
   //
   //uint64_t startTime = ofGetElapsedTimeMicros();
   //Cost total_cost = hungarian_algorithm::solve<Cost>(cost_matrix, cost_matrix.rows(), cost_matrix.cols(), assignment);
   //uint64_t elapsedTime = ofGetElapsedTimeMicros() - startTime;

   //cout << "Assignments are: " << endl;
   //for (int i = 0; i < size; i++) {
   //    cout << "\t" << i << "\t:\t" << assignment[i] << endl;
   //}
   ////cout << "Assignments are: " << assignment[0] << "\t" << assignment[1] << endl;
   //cout << "Total cost = " << total_cost << endl;
   //cout << "Took " << float(elapsedTime) / 1000.0 << " ms" << endl;





   //typedef float Cost;

   //// Rows represents the new sample; columns represent the old sample
   //Matrix2D<Cost> cost_matrix(2, 3);
   //cost_matrix(0, 0) = -1;
   //cost_matrix(0, 1) = 0;
   //cost_matrix(0, 2) = 0;
   //cost_matrix(1, 0) = 0;
   //cost_matrix(1, 1) = -1;
   //cost_matrix(1, 2) = 0;
   //cout << cost_matrix << endl;

   //vector<size_t> assignment(cost_matrix.rows());

   //uint64_t startTime = ofGetElapsedTimeMicros();
   //Cost total_cost = hungarian_algorithm::solve<Cost>(cost_matrix, cost_matrix.rows(), cost_matrix.cols(), assignment);
   //uint64_t elapsedTime = ofGetElapsedTimeMicros() - startTime;

   //// For each new sample (row), what's the matching old sample (column)?
   //cout << "Assignments are: " << endl;
   //for (int i = 0; i < cost_matrix.rows(); i++) {
   //    cout << "\t" << i << "\t:\t" << assignment[i] << endl;
   //}
   ////cout << "Assignments are: " << assignment[0] << "\t" << assignment[1] << endl;
   //cout << "Total cost = " << total_cost << endl;
   //cout << "Took " << float(elapsedTime) / 1000.0 << " ms" << endl;




    typedef float Cost;

    // Rows represents the new sample; columns represent the old sample
    Matrix2D<Cost> cost_matrix(3, 2);
    cost_matrix(0, 0) = 1;
    cost_matrix(0, 1) = 0;
    cost_matrix(1, 0) = 0;
    cost_matrix(1, 1) = 1;
    cost_matrix(2, 0) = 1;
    cost_matrix(2, 1) = 1;
    cout << cost_matrix << endl;

    vector<size_t> assignment(cost_matrix.rows());

    uint64_t startTime = ofGetElapsedTimeMicros();
    Cost total_cost = hungarian_algorithm::solve<Cost>(cost_matrix, cost_matrix.rows(), cost_matrix.cols(), assignment);
    uint64_t elapsedTime = ofGetElapsedTimeMicros() - startTime;

    // For each new sample (row), what's the matching old sample (column)?
    cout << "Assignments are: " << endl;
    for (int i = 0; i < cost_matrix.rows(); i++) {
        cout << "\t" << i << "\t:\t" << assignment[i] << endl;
    }
    //cout << "Assignments are: " << assignment[0] << "\t" << assignment[1] << endl;
    cout << "Total cost = " << total_cost << endl;
    cout << "Took " << float(elapsedTime) / 1000.0 << " ms" << endl;

}

//--------------------------------------------------------------
void ofApp::update(){
    
}

//--------------------------------------------------------------
void ofApp::draw(){

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
