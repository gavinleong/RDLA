#include <boost/function_output_iterator.hpp>
#include <boost/geometry/geometry.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <ostream>
#include <numeric>
#include <algorithm>

// using cxxopts for CLI argument parsing
#include "packages/cxxopts.hpp"

// using prakhar1989's ProgressBar for CLI progress indicator
#include "packages/ProgressBar.hpp"

// using tinyobjloader for OBJ file parsing
#define TINYOBJLOADER_IMPLEMENTATION
#define EPSILON 0.00000001
#include "packages/tiny_obj_loader.h"

// using Eigen for matrix operations
#include <Eigen/Dense>

using namespace std;

// full 3D model to inhibit growth
//string fullModelFilename = "1024_596_fixed_box.obj";
string fullModelFilename = "F595_flat_s53_sfm_processed_box.obj";
tinyobj::attrib_t fullModelAttrib;
vector<tinyobj::shape_t> fullModelShapes;
vector<tinyobj::material_t> fullModelMaterials;

// 3D model of selected faces to seed growth
//string seedModelFilename = "1024_596_fixed.obj";
string seedModelFilename = "F595_flat_s53_sfm_processed.obj";
tinyobj::attrib_t seedModelAttrib;
vector<tinyobj::shape_t> seedModelShapes;
vector<tinyobj::material_t> seedModelMaterials;

// number of particles



const int DefaultNumberOfParticles = 10000;
int numParticles = DefaultNumberOfParticles;

// height of RDLA
const double DefaultRdlaHeight = 0.010; // in millimetres
double RdlaHeight = DefaultRdlaHeight;

// progress bar
ProgressBar progressBar(numParticles, 40);

// output file
string filename = "Cone_parentcorrected";
ofstream file;


// percentage at which attractionDistance begins to increase 
double attractionBoost = 0.95;

// id number for the m_IndexSmall final index for searching over
int id2 = 0;

// for running the IndexSmall loop for collecting the slices not yet full 
// once the loop has run once, this will turn to True.
bool IndexSmallDone = false;

// earlyStop prevents the DLA from running until forever before we get the full DLA histogram effect.
// From testing with F596, we find that the slow-down starts around 95% completion, but just in case, 
// as smaller meshes might slow down earlier due to not having enough points to attach to, we'll use 
// 93% stop as a precaution to prevent getting stuck.
double earlyStop = 1.0;

// for the findIndex() function, if it becomes 1, I have set a trigger pause condition during debug
int index_error = 0;

// knn_slice sets the number of nearest neighbouring points to a query point for the last 1- fastSearchModePerc
// percentage the RDLA. This may speed up the last bit of RDLA as it actively searches for viable points to 
// stick to.
int knn_slice = 25;

// fast search mode turns on towards the end of the RDLA run, so that for the second-to-last stage the points do not
// take a long time to attach to a viable parent.
double fastSearchModePerc = 0.80;

// fast search mode turns on for the last few iterations towards the end of the RDLA run, so that for the
// last stage the points do not take a long time to attach to a viable parent (alternative to regular DLA mode).
double smallIndexPerc = 0.80;

// DLA mode turns as the last stage of the RDLA run, so that the last points do not ground 
// the RDLA to a snail-like crawl for attaching to a viable parent.
double DLAModePerc = 0.95;

// counter for the number of points currently in each slice
unsigned int ptsPerSliceCounter[31] = { };

// interval (in iterations) that point data is outputted
int lastOutputIteration = 0;
int interval = numParticles;

// number of dimensions (must be 2 or 3)
const int D = 3;

// default parameters (documented below)
// default area of the seed mesh (needs to be defined externally with a .csv file 
// calculated from the compute Geometric Measures from MeshLab
const double DefaultSeedArea = 0.031138701;

double SeedArea = DefaultSeedArea;

const double DefaultParticleSpacing = 0.0012;
const double DefaultAttractionDistance = 0.003;
const double DefaultMinMoveDistance = 0.01;
const int DefaultStubbornness = 1;
const double DefaultStickiness = 1.0;
const double DefaultBoundingRadius = 0.01;
// Angle in degrees for the cone width (defines the angle between the (0,0,1) and the max
// angle: DefaultConeAngleDegree
const double DefaultConeAngleDegree = 45.0;
// DefaultConeReductionRate is between 0 to 1, where 0 is slow cone reduction rate and 1 is
// fast cone reduction rate. A slow cone reduction means the DLA structure can bend and curl
// for more generations before it “decides” on a direction. Fast means the DLA structure 
// chooses direction in fewer generations, resulting in perhaps more angular DLA, less bending.
const double DefaultConeReductionRate = 0.2;
// DistanceDiffusionStarts is a double between 0.0 and 1.0, where the smaller the number,
// the less likely the DLA will start earlier on in the structure, and the larger, the
// more likely DLA will start early in the structure.

// choose value between 0 and 0.0375 (which is in metres), for the distance between the rock
// surface and a point above the rock surface where you want the mean value of attachment mode,
// i.e., value below which there is 50% chance of the point sticking mode being cone-based and
// value above which there is 50% chance of RDLA-based attachment. 
const double DefaultDistanceDiffusionStarts = 0.01;


// define the thickness slices, matching the real lichen laser histogram with 31 bins
double thicknessRef[32] = {0.0, 0.00125, 0.00250, 0.00375, 0.00500, 0.00625, 0.00750, 0.00875, 0.01000, 0.01125,
            0.01250, 0.01375, 0.01500, 0.01625, 0.01750, 0.01875, 0.02000, 0.02125, 0.02250, 0.02375, 0.02500,
            0.02625, 0.02750, 0.02875, 0.03000, 0.03125, 0.03250, 0.03375, 0.03500, 0.03625, 0.03750, 1.0};

// Constant to add
double sheathThickness = 0.00085; // for adjusting for the sheathing via distance mapping in Python

// define the probability per bin, which is proportionate to the height of the histogram
double thicknessProbability[31] = {1914, 2201, 2176, 2227, 2294, 2461, 2500, 2271, 2221, 1957,
            1597, 1347, 1275, 1086, 866, 764, 546, 452, 347, 276, 179, 122, 91, 61, 29, 14, 20, 11, 5,
            4, 0};

//double thicknessProbability[31] = {0.06665195, 0.057386, 0.06172908, 0.08468703, 0.08861486, 0.09967235,
//0.09618321, 0.08250105, 0.07449582, 0.0586262, 0.03999661, 0.03388108,
//0.03397434, 0.03089567, 0.02181267, 0.02093657, 0.01360918, 0.01092247,
//0.00798724, 0.00640057, 0.00294448, 0.00212911, 0.00134974, 0.00104483,
//0.00025328, 0.00032925, 0.00108742, 0.00058987, 0.00028152, 0.00019161, 0};

// total probability to use as a normalisation term for calculating the maximum number of points per slice
double totalProbability = accumulate(begin(thicknessProbability), end(thicknessProbability), 0.0, plus<double>());
// number of points allowed for each slice based on probability and total number of points in simulation
vector<double> ptsPerSlice(31);



// scale value for the logistic function, for adapting the function to a range of input values
const double logisticFuncScale = 0.002;

// vector for storing all the distance value of the attached points to the mesh. For importing into
// python code for calculating histogram
vector<double> distanceVector(numParticles);

// boost is used for its spatial index
using BoostPoint = boost::geometry::model::point<double, D, boost::geometry::cs::cartesian>;
using IndexValue = pair<BoostPoint, int>;
using Index = boost::geometry::index::rtree<IndexValue, boost::geometry::index::linear<4>>;

// VARIABLE FOR PRINTING THE PARENT PARTICLE LOCATION
//int parentRow = 0;
//int parentRowIndex[DefaultNumberOfParticles];
//int parentRowIndexCounter = 0;
//double parent_vector[DefaultNumberOfParticles][3];

// VARIABLE FOR PRINTING THE FACE PARTICLE LOCATION
//int faceRow = 0;
//int faceRowIndex[DefaultNumberOfParticles];
//int faceRowIndexCounter = 0;
//double face_vector[DefaultNumberOfParticles][3];

// DEFINE PI    
constexpr double pi = 3.14159265358979323846;


// Vector represents a point or a vector
class Vector {
public:
    Vector() :
        m_X(0), m_Y(0), m_Z(0) {}

    Vector(double x, double y) :
        m_X(x), m_Y(y), m_Z(0) {}

    Vector(double x, double y, double z) :
        m_X(x), m_Y(y), m_Z(z) {}

    double X() const {
        return m_X;
    }

    double Y() const {
        return m_Y;
    }

    double Z() const {
        return m_Z;
    }

    BoostPoint ToBoost() const {
        return BoostPoint(m_X, m_Y, m_Z);
    }

    double Length() const {
        return sqrt(m_X * m_X + m_Y * m_Y + m_Z * m_Z);
    }

    double LengthSquared() const {
        return m_X * m_X + m_Y * m_Y + m_Z * m_Z;
    }

    double Distance(const Vector& v) const {
        const double dx = m_X - v.m_X;
        const double dy = m_Y - v.m_Y;
        const double dz = m_Z - v.m_Z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    Vector Normalized() const {
        const double m = 1 / Length();
        return Vector(m_X * m, m_Y * m, m_Z * m);
    }

    Vector GetCrossed(const Vector& v) const {
        return Vector(m_Y * v.m_Z - m_Z * v.m_Y,
            m_Z * v.m_X - m_X * v.m_Z,
            m_X * v.m_Y - m_Y * v.m_X);
    }

    double GetDot(const Vector& v) const {
        return m_X * v.m_X +
            m_Y * v.m_Y +
            m_Z * v.m_Z;
    }

    Vector operator+(const Vector& v) const {
        return Vector(m_X + v.m_X, m_Y + v.m_Y, m_Z + v.m_Z);
    }

    Vector operator-(const Vector& v) const {
        return Vector(m_X - v.m_X, m_Y - v.m_Y, m_Z - v.m_Z);
    }

    Vector operator*(const double a) const {
        return Vector(m_X * a, m_Y * a, m_Z * a);
    }

    Vector& operator+=(const Vector& v) {
        m_X += v.m_X; m_Y += v.m_Y; m_Z += v.m_Z;
        return *this;
    }

    bool operator==(const Vector& v) const {
        return m_X == v.m_X && m_Y == v.m_Y && m_Z == v.m_Z;
    }

private:
    double m_X;
    double m_Y;
    double m_Z;
};


// vector of coordinates of the attached points, put into a vector the same size as ptsPerSlice. Each of 
// the bins of pointsSlice contains the coordinates of the points that the slice contains.
vector<vector<Vector>> pointsSlice(31, vector<Vector>(0));

// Lerp linearly interpolates from a to b by distance.
Vector Lerp(const Vector& a, const Vector& b, const double d) {
    return a + (b - a).Normalized() * d;
}

// Random returns a uniformly distributed random number between lo and hi
double Random(const double lo = 0, const double hi = 1) {
    static thread_local mt19937 gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> dist(lo, hi);
    return dist(gen);
}





///////////////////////////////////////////////////////////////////////////////////////////////////
// Gavin custom:
// Returns a uniformly distributed random integer between lo and hi (inclusive)
int RandomInteger(const int lo, const int hi) {
    static thread_local mt19937 gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> dist(lo, hi);
    return dist(gen);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// function for getting current date and time, for printing to file name 
const string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    errno_t timenow = localtime_s(&tstruct, &now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%H-%M-%S", &tstruct);
    return buf;
}
///////////////////////////////////////////////////////////////////////////////////////////////////



// RandomInUnitSphere returns a random, uniformly distributed point inside the
// unit sphere (radius = 1)
Vector RandomInUnitSphere() {
    while (true) {
        const Vector p = Vector(
            Random(-1, 1),
            Random(-1, 1),
            D == 2 ? 0 : Random(-1, 1));
        if (p.LengthSquared() < 1) {
            return p;
        }
    }
}


//custom code below
////////////////////////////////////////////////////////////////////////////////////////////////////////////



// MaxMeshBound returns a vector for the maximum x, y and z of a mesh
Vector MaxMeshBound() {

    vector<Vector> vertices;
    vertices.push_back(Vector(0.0, 0.0, 0.0));

    size_t index_offset = 0;

    // go over all the faces in the shape
    for (size_t f = 0; f < fullModelShapes[0].mesh.num_face_vertices.size(); f++) {
        unsigned int fv = fullModelShapes[0].mesh.num_face_vertices[f];

        // collect the vertices
        for (size_t v = 0; v < fv; v++) {
            tinyobj::index_t idx = fullModelShapes[0].mesh.indices[v + index_offset];
            tinyobj::real_t vx = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 0];
            tinyobj::real_t vy = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 1];
            tinyobj::real_t vz = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 2];

            //cout << endl << static_cast<double>(vertices[0].X()) << endl;
            //cout << endl << static_cast<double>(vertices[0].Y()) << endl;
            //cout << endl << static_cast<double>(vertices[0].Z()) << endl;
            vertices[0] = Vector(max(vertices[0].X(), static_cast<double>(vx)),
                max(vertices[0].Y(), static_cast<double>(vy)),
                max(vertices[0].Z(), static_cast<double>(vz)));
        }
        index_offset += fv;
    }
    return vertices[0];
}


Vector MinMeshBound() {

    vector<Vector> vertices;
    vertices.push_back(Vector(0.0, 0.0, 0.0));

    size_t index_offset = 0;

    // go over all the faces in the shape
    for (size_t f = 0; f < fullModelShapes[0].mesh.num_face_vertices.size(); f++) {
        unsigned int fv = fullModelShapes[0].mesh.num_face_vertices[f];

        // collect the vertices
        for (size_t v = 0; v < fv; v++) {
            tinyobj::index_t idx = fullModelShapes[0].mesh.indices[v + index_offset];
            tinyobj::real_t vx = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 0];
            tinyobj::real_t vy = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 1];
            tinyobj::real_t vz = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 2];

            //cout << endl << static_cast<double>(vertices[0].X()) << endl;
            //cout << endl << static_cast<double>(vertices[0].Y()) << endl;
            //cout << endl << static_cast<double>(vertices[0].Z()) << endl;
            vertices[0] = Vector(min(vertices[0].X(), static_cast<double>(vx)),
                min(vertices[0].Y(), static_cast<double>(vy)),
                min(vertices[0].Z(), static_cast<double>(vz)));
        }
        index_offset += fv;
    }
    return vertices[0];
}


// the following function is custom code that calculates the area of the seed mesh by
// calculating the area of each of its triangles via the cross product and adding them
// up. Will be used to calculate gamma and adapt the height of RDLA 


// Calculate the area of each face of a mesh and return the total area
double CalculateTotalFaceArea() {

    double totalArea = 0.0;
    ofstream file;
    file.open("C:\\Users\\gavin\\OneDrive - University of Brighton\\All_UCL_onedrive\\ALL_CARVINGS_4_MESHNET\\NON-CARVINGS\\check595.csv");


    // go over all the faces in the shape
    for (size_t f = 0; f < seedModelShapes[0].mesh.num_face_vertices.size(); f++) {
        unsigned int fv = seedModelShapes[0].mesh.num_face_vertices[f];

        // Check if the face is a triangle (3 vertices) to calculate its area
        if (fv == 3) {
            tinyobj::index_t idx1 = seedModelShapes[0].mesh.indices[f * 3 + 0];
            tinyobj::index_t idx2 = seedModelShapes[0].mesh.indices[f * 3 + 1];
            tinyobj::index_t idx3 = seedModelShapes[0].mesh.indices[f * 3 + 2];

            tinyobj::real_t vx1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 0];
            tinyobj::real_t vy1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 1];
            tinyobj::real_t vz1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 2];

            tinyobj::real_t vx2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 0];
            tinyobj::real_t vy2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 1];
            tinyobj::real_t vz2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 2];

            tinyobj::real_t vx3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 0];
            tinyobj::real_t vy3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 1];
            tinyobj::real_t vz3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 2];

            // Calculate the vectors formed by the vertices of the face
            Vector v1(vx1, vy1, vz1);
            Vector v2(vx2, vy2, vz2);
            Vector v3(vx3, vy3, vz3);

            // Calculate the cross product to get the face normal and its magnitude to get the area
            Vector crossProduct = (v2 - v1).GetCrossed(v3 - v1);
            double area = 0.5 * crossProduct.Length();
            //cout << area << endl;
            //cout << area << endl;


            
            if (file.is_open()) {
                file << area << endl;
                cout << "open" << endl;
            }


            // Add the area of this face to the total area
            totalArea += area;

        }
        // For faces with more than 3 vertices (non-triangular faces), you need more complex methods to calculate the area.
        // This example assumes you only have triangular faces.
    }
    file.close();
    return totalArea;
}




/*
// Calculate the area of each face of a mesh and return the total area
double CalculateTotalFaceAreaPrint() {

    double totalArea = 0.0;

    ofstream file;

    file.open("data/reprintobj_3.csv");

    if (file.is_open()) {
        cout << endl << "file opened";
    }
    // go over all the faces in the shape
    for (size_t f = 0; f < seedModelShapes[0].mesh.num_face_vertices.size(); f++) {
        unsigned int fv = seedModelShapes[0].mesh.num_face_vertices[f];

        // Check if the face is a triangle (3 vertices) to calculate its area
        if (fv == 3) {
            tinyobj::index_t idx1 = seedModelShapes[0].mesh.indices[f * 3 + 0];
            tinyobj::index_t idx2 = seedModelShapes[0].mesh.indices[f * 3 + 1];
            tinyobj::index_t idx3 = seedModelShapes[0].mesh.indices[f * 3 + 2];

            cout << "f: " << f << endl;
            cout << "id1: " << idx1.vertex_index << endl;
            cout << "id2: " << idx2.vertex_index << endl;
            cout << "id3: " << idx3.vertex_index << endl;

            tinyobj::real_t vx1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 0];
            tinyobj::real_t vy1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 1];
            tinyobj::real_t vz1 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx1.vertex_index + 2];

            tinyobj::real_t vx2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 0];
            tinyobj::real_t vy2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 1];
            tinyobj::real_t vz2 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx2.vertex_index + 2];

            tinyobj::real_t vx3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 0];
            tinyobj::real_t vy3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 1];
            tinyobj::real_t vz3 = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx3.vertex_index + 2];

            // Calculate the vectors formed by the vertices of the face
            Vector v1(vx1, vy1, vz1);
            Vector v2(vx2, vy2, vz2);
            Vector v3(vx3, vy3, vz3);

            if (f == 0) {
                file << v1.X() << "," << v1.Y() << "," << v1.Z() << endl;
                file << v2.X() << "," << v2.Y() << "," << v2.Z() << endl;
                file << v3.X() << "," << v3.Y() << "," << v3.Z() << endl;
            }

            // Calculate the cross product to get the face normal and its magnitude to get the area
            Vector crossProduct = (v2 - v1).GetCrossed(v3 - v1);
            double area = 0.5 * crossProduct.Length();
            //cout << area << endl;

            // Add the area of this face to the total area
            totalArea += area;
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////
    }

    // For faces with more than 3 vertices (non-triangular faces), you need more complex methods to calculate the area.
    // This example assumes you only have triangular faces.
    file.close();
    return totalArea;
}


*/


//Vector CentreOfMass() {
//
//    vector<Vector> vertices;
//    vertices.push_back(Vector(0.0, 0.0, 0.0));
//
//    size_t index_offset = 0;
//
//    // go over all the faces in the shape
//    for (size_t f = 0; f < fullModelShapes[0].mesh.num_face_vertices.size(); f++) {
//        unsigned int fv = fullModelShapes[0].mesh.num_face_vertices[f];
//
//        // collect the vertices
//        for (size_t v = 0; v < fv; v++) {
//            tinyobj::index_t idx = fullModelShapes[0].mesh.indices[v + index_offset];
//            tinyobj::real_t vx = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 0];
//            tinyobj::real_t vy = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 1];
//            tinyobj::real_t vz = fullModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 2];
//
//            //cout << endl << static_cast<double>(vertices[0].X()) << endl;
//            //cout << endl << static_cast<double>(vertices[0].Y()) << endl;
//            //cout << endl << static_cast<double>(vertices[0].Z()) << endl;
//            vertices[0] = Vector(min(vertices[0].X(), static_cast<double>(vx)),
//                min(vertices[0].Y(), static_cast<double>(vy)),
//                min(vertices[0].Z(), static_cast<double>(vz)));
//        }
//        index_offset += fv;
//    }
//    return vertices[0];
//}



// INCOMPLETE CODE, MAKE IT CHOOSE POINT USING MAX AND MIN MESHBOUND() FUNC.

// RandomInMesh returns a random, uniformly distributed point inside the
// mesh bounds
Vector RandomInMesh() {
    while (true) {
        vector<Vector> minBound;
        vector<Vector> maxBound;
        minBound.push_back(MinMeshBound());
        maxBound.push_back(MaxMeshBound());
        const Vector p = Vector(
            Random(static_cast<double>(minBound[0].X()), static_cast<double>(maxBound[0].X())),
            Random(static_cast<double>(minBound[0].Y()), static_cast<double>(maxBound[0].Y())),
            D == 2 ? 0 : Random(static_cast<double>(minBound[0].Z()), static_cast<double>(maxBound[0].Z())));
        if (p.LengthSquared() < max(minBound[0].LengthSquared(), maxBound[0].LengthSquared())) {
            return p;
        }
    }
}


// calculates the gamma, which is the conversion scalar to get from RDLA height and area to 
// number of points n.
/*
void GammaCalculate() {

    double meshArea = CalculateTotalFaceArea();

    double gamma = RdlaHeight / (numParticles * meshArea);

    cout << gamma;
}
*/

// custom code: compares one variable, say x, to an array of 31 values, say D,
// and finds the index of D where x is larger than the value of D at that index, 
// but smaller than the value at that index + 1
int findIndex(double x, double D[], int size) {
    for (int i = 0; i < size - 1; ++i) {
        if (x >= D[i] && x < D[i + 1]) {
            //if (D[i + 1] == 1.0) {
                //cout << "point is between " << D[i] << " and " << D[i + 1] << endl;
            //}
            return i;
        }
    }
    index_error = 1;
    cout << "findIndex not working" << endl;

    return -1; // Return -1 if x doesn't fit within the range of the array

}



// custom code: compares one variable, say x, to an array of 31 values, say D,
// and finds the index of D where x is larger than the value of D at that index, 
// but smaller than the value at that index + 1
unsigned int findIndex2(double x, double D[], int size) {
    for (int i = 0; i < size - 1; ++i) {
        if (x >= D[i] && x < D[i + 1]) {
            if (D[i + 1] == 1.0) {
                cout << D[i + 1] << endl;
            }
            return i;
        }
    }
    return -1; // Return -1 if x doesn't fit within the range of the array
}


// returns string with defined decimal precision for saving parameters to file names 
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}



//custom code above
////////////////////////////////////////////////////////////////////////////////////////////////////////////




// IsIntersectingFace checks for intersection of a ray (O) with random direction (D) and a triangle face defined by vertices V1, V2, and V3
// Adapted from Amnon Owed's ofxPointInMesh::triangleIntersection: https://github.com/AmnonOwed/ofxPointInMesh
// Uses Möller–Trumbore: http://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool IsIntersectingFace(const Vector& V1, const Vector& V2, const Vector& V3, const Vector& O, const Vector& D, Vector& R) {
    Vector e1, e2; // Edge1, Edge2
    Vector P, Q, T;
    float det, inv_det, u, v;
    float t;

    // Find vectors for two edges sharing V1
    e1 = V2 - V1;
    e2 = V3 - V1;

    // Begin calculating determinant - also used to calculate u parameter
    P = D.GetCrossed(e2);

    // if determinant is near zero, ray lies in plane of triangle
    det = e1.GetDot(P);

    // NOT CULLING
    if (det > -EPSILON && det < EPSILON) {
        return false;
    }

    inv_det = 1.f / det;

    // calculate distance from V1 to ray origin
    T = O - V1;

    // calculate u parameter and test bound
    u = T.GetDot(P) * inv_det;

    // the intersection lies outside of the triangle
    if (u < 0.f || u > 1.f) {
        return false;
    }

    // prepare to test v parameter
    Q = T.GetCrossed(e1);

    // calculate V parameter and test bound
    v = D.GetDot(Q) * inv_det;

    // the intersection lies outside of the triangle
    if (v < 0.f || u + v  > 1.f) {
        return false;
    }

    t = e2.GetDot(Q) * inv_det;

    if (t > EPSILON) { // ray intersection
        R = O + D * t; // store intersection point
        return true;
    }

    // no hit, no win
    return false;
}

// IsInsideMesh determines if a given point is inside of the base model mesh
// Adapted from Amnon Owed's ofxPointInMesh::isInside: https://github.com/AmnonOwed/ofxPointInMesh
bool IsInsideMesh(Vector& p) {
    // if no mesh was provided, skip this test
    if (fullModelFilename.empty()) {
        return false;
    }

    Vector foundIntersection; // variable to store a single found intersection
    vector<Vector> results;  // vector to store all found intersections
    Vector randomDirection = Vector(0.1, 0.2, 0.3); // a random direction

    // go over all the shapes in the mesh
    for (size_t s = 0; s < fullModelShapes.size(); s++) {
        size_t index_offset = 0;

        // go over all the faces in the shape
        for (size_t f = 0; f < fullModelShapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = fullModelShapes[s].mesh.num_face_vertices[f];
            vector<Vector> vertices;

            // collect the vertices
            for (size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = fullModelShapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = fullModelAttrib.vertices[3 * idx.vertex_index + 0];
                tinyobj::real_t vy = fullModelAttrib.vertices[3 * idx.vertex_index + 1];
                tinyobj::real_t vz = fullModelAttrib.vertices[3 * idx.vertex_index + 2];
                vertices.push_back(Vector(vx, vy, vz));
            }

            index_offset += fv;

            // do a triangle-ray intersection on each face in the mesh
            // store the intersection (if any) in the variable foundIntersection
            if (IsIntersectingFace(vertices[0], vertices[1], vertices[2], p, randomDirection, foundIntersection)) {
                // store all found intersections
                results.push_back(foundIntersection);
            }
        }
    }

    // handle multiple mesh intersections at the same point (by removing duplicates)
    vector<Vector> unique_results;
    unique_copy(results.begin(), results.end(), back_inserter(unique_results));

    // // determine if the point is inside or outside the mesh, based on the number of unique intersections
    if (unique_results.size() % 2 == 1) {
        return true;
    }
    else {
        return false;
    }
}

// GetDot2 calculates the dot product of a vector with itself
double GetDot2(const Vector& p) {
    return p.GetDot(p);
}

// GetSign returns -1 for negative numbers, 1 for positive numbers, and 0 for 0
int GetSign(const double& v) {
    int result = 0;
    if (v < 0) {
        result = -1;
    }
    else if (v > 0) {
        result = 1;
    }
    return result;
}

// clamp ensures that a value (x) is within the range defined by [lower] and [upper]
double clamp(double x, double lower, double upper) {
    return min(upper, max(x, lower));
}

// DistanceToTriangle calculates the shortest distance between a point (p) and a 3D triangle defined by vertices v1, v2, and v3
// Adapted from: https://iquilezles.org/www/articles/triangledistance/triangledistance.htm
double DistanceToTriangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& p) {
    Vector v21 = v2 - v1;   
    Vector v32 = v3 - v2;
    Vector v13 = v1 - v3;
    Vector p1 = p - v1;
    Vector p2 = p - v2;
    Vector p3 = p - v3;
    Vector nor = v21.GetCrossed(v13);

    return sqrt(
        // inside/outside test
        (GetSign(p1.GetDot(v21.GetCrossed(nor))) +
            GetSign(p2.GetDot(v32.GetCrossed(nor))) +
            GetSign(p3.GetDot(v13.GetCrossed(nor))) < 2.0)
        ?
        // 3 edges
        min(
            min(
                GetDot2(v21 * clamp(p1.GetDot(v21) / GetDot2(v21), 0.0, 1.0) - p1),
                GetDot2(v32 * clamp(p2.GetDot(v32) / GetDot2(v32), 0.0, 1.0) - p2)
            ),
            GetDot2(v13 * clamp(p3.GetDot(v13) / GetDot2(v13), 0.0, 1.0) - p3)
        )
        :
        // 1 face
        p1.GetDot(nor) * p1.GetDot(nor) / GetDot2(nor)
    );
}


// Model holds all of the particles and defines their behavior.
class Model {
public:
    Model() :
        m_ParticleSpacing(DefaultParticleSpacing),
        m_AttractionDistance(DefaultAttractionDistance),
        m_MinMoveDistance(DefaultMinMoveDistance),
        m_Stubbornness(DefaultStubbornness),
        m_Stickiness(DefaultStickiness),
        m_BoundingRadius(DefaultBoundingRadius),
        m_ConeAngleDegree(DefaultConeAngleDegree),
        m_ConeReductionRate(DefaultConeReductionRate),
        m_DistanceDiffusionStarts(DefaultDistanceDiffusionStarts) {}
        //m_SeedArea(DefaultSeedArea) {}

    void SetParticleSpacing(const double a) {
        m_ParticleSpacing = a;
    }

    void SetAttractionDistance(const double a) {
        m_AttractionDistance = a;
    }

    void SetMinMoveDistance(const double a) {
        m_MinMoveDistance = a;
    }

    void SetStubbornness(const int a) {
        m_Stubbornness = a;
    }

    void SetStickiness(const double a) {
        m_Stickiness = a;
    }

    void SetBoundingRadius(const double a) {
        m_BoundingRadius = a;
    }

    void SetConeReductionRate(const double a) {
        m_ConeReductionRate = a;
    }

    void SetConeAngleDegree(const double a) {
        m_ConeAngleDegree = a;
    }

    void SetDistanceDiffusionStarts(const double a) {
        m_DistanceDiffusionStarts = a;
    }

    //void SetSeedArea(const double a) {
    //    m_SeedArea = a;
    //}

    // returns a direction at which a parent point can attach a new point to it
    // in a cone of directions. See links below for source (converted from MATLAB code)
    // https://stackoverflow.com/questions/38997302/create-random-unit-vector-inside-a-defined-conical-region
    // https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone/205589#205589
    Vector RandomDirectionInCone(const int parent) {

        //if (check == 0)
        //    return Vector();

        //variables for randomised cone vertex
        double xRot;
        double yRot;
        Eigen::Vector3d r;
        double zRot;
        double phi;

        //// NEED TO ADD A CONDITION FOR m_GenParent above 1 to move to defining triangleNormal (change its name)
        //// to the previous direction of RandomDirectionInCone return
        //if (m_GenerationParent > 1) {
        //    Vector triangleNormal = 
        //}
        //else {
        Vector structureDirection = m_StructureDirection[parent];
        //}
        //Vector structureDirection = m_StructureDirection[m_StructureDirection.size() - 1];

        //variables for translating rotation to structureDirection axes
        Eigen::Matrix3d R = Eigen::Array33d::Zero();

        /////////////////////////////////////////////////////////////
        //identity matrix 
        //vector<vector<int>> eye;

        //// Initialize the matrix as a n x n array of 0.
        //eye = vector<vector<int>>(3, vector<int>(3, 0));

        //// Set the diagonal to be 1s
        //for (unsigned int t = 0; t < 3; t++)
        //    eye[t][t] = 1;
        /////////////////////////////////////////////////////////////
        Eigen::Matrix3d eye;
        eye << 1, 0, 0,
            0, 1, 0,
            0, 0, 1;


        // random point on spherical cap in direction of vector (0,0,1)
        // [1] See https ://math.stackexchange.com/a/205589/81266

        // exponential decay equation used as method of reducing the m_ConeAngleDegree
        // m_ConeReductionRate is between 0 and 1: the larger the number, the fewer generations it takes to 
        // reduce m_ConeAngleDegree to small values and so the thallus "straightens".
        // m_ConeAngleDegree is the starting angle value.
        // since we set the first m_GenerationParent, i.e. point on the mesh face, as = 1, we must subtract 1 from 
        // the parent to let the first coneAngle value = m_ConeAngleDegree.
        // I am removing the 1 from m_GenerationParent[parent - 1] because I think the -1 is not supposed to be there. 
        // parent can be 0, so that means it can be -1, but if it's -1, then it goes out of the index range of m_GenerationParent
        // since that array cannot go to -1.
        double coneAngle = (m_ConeAngleDegree * exp(-m_ConeReductionRate * m_GenerationParent[parent])) * pi / 180;

        zRot = Random() * (1 - cos(coneAngle)) + cos(coneAngle);
        phi = Random() * 2 * pi;
        xRot = sqrt(1 - (zRot * zRot)) * cos(phi);
        yRot = sqrt(1 - (zRot * zRot)) * sin(phi);


        //do {
        //    xRot = 2.0 * Random() - 1.0;
        //    yRot = 2.0 * Random() - 1.0;
        //} while (xRot * xRot + yRot * yRot > 1.0);
        //// (x,y) is now randomly distributed within the unit circle

        ////theta is maximum angle from the z axis for the cone width
        //r = tan(m_Theta) * m_ConeReductionRate;
        //Vector randomInCone(r * xRot, r * yRot, 1);


        // converted from matlab to c++, but need to check it is correct
        //Find the rotation axis `u` and rotation angle `rot`
        Vector u = (Vector(0, 0, 1).GetCrossed(structureDirection)).Normalized();
        Eigen::Vector3d v(u.X(), u.Y(), u.Z());
        double rot = acos(structureDirection.GetDot(Vector(0, 0, 1)));

        //Convert rotation axis and angle to 3x3 rotation matrix
        //See https ://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
        //double crossMatrix[3][3] = { {0, -u.Z(), u.Y()}, {u.Z(), 0, -u.X()}, {-u.Y(), u.X(), 0} };
        // aka the Rodrigues' rotation formula (which requires the cross-product matrix): 
        // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        Eigen::Matrix3d crossMatrix;
        crossMatrix << 0, -u.Z(), u.Y(),
            u.Z(), 0, -u.X(),
            -u.Y(), u.X(), 0;
        R = cos(rot) * eye + sin(rot) * crossMatrix + (1 - cos(rot)) * (v * v.adjoint());

        Eigen::Vector3d xyz(xRot, yRot, zRot);

        //Rotate[x; y; z] from north pole to `coneDir`.
        r = R * xyz;
        //double r_access = double(r[0]);
        //vector<double> v2;
        //v2.resize(r.size());
        //Eigen::Vector3d::Map(&v2[0], r.size()) = r;

        Vector direction = Vector(r[0], r[1], r[2]);

        // stores the SeedstructureDirection of the parent particle for particle attaching to parent
        m_StructureDirection.push_back(direction.Normalized());
        //cout  << m_StructureDirection[m_StructureDirection.size() - 1].X() << ", "
        //   << m_StructureDirection[m_StructureDirection.size() - 1].Y() << ", "
        //  << m_StructureDirection[m_StructureDirection.size() - 1].Z() << endl;
        /////////////////////////////////////////////////////////////////////////

        return direction;
    }




    // Add adds a new particle with the specified parent particle, this is used when a 
    // particle attaches to a parent particle
    void Add(const Vector& p, const int parent = -1) {
        const int id = m_Points.size();
        m_Index.insert(make_pair(p.ToBoost(), id));
        m_Points.push_back(p);
        m_Parents.push_back(parent);
        m_JoinAttempts.push_back(0);
        m_BoundingRadius = max(m_BoundingRadius, p.Length() + m_AttractionDistance);

        /////////////////////////////////////////////////////////////////////////
        // for particle attaching to a parent, increment the m_GenerationParent of parent to get
        // m_GenerationParent of the current particle.
        m_GenerationParent.push_back(m_GenerationParent[parent]++);
        //cout << endl << "generation: " << m_GenerationParent[m_GenerationParent.size() - 1];

        // THIS m_StructureDirection DEFINITION IS MOVED TO THE RandomDirectionInCone FUNCTION
        // stores the StructureDirection of the parent particle for particle attaching to parent
        //m_StructureDirection.push_back(m_StructureDirection[parent]);
        //cout << endl << m_StructureDirection[m_StructureDirection.size() - 1].X() << ", "
        //    << m_StructureDirection[m_StructureDirection.size() - 1].Y() << ", "
        //    << m_StructureDirection[m_StructureDirection.size() - 1].Z();
        /////////////////////////////////////////////////////////////////////////

        // update and display progress bar
        ++progressBar;
        progressBar.display();

        // wrap up the progress bar when last particle is placed
        if (id == numParticles) {
            progressBar.done();
        }
    }

    // Function overloading of Add()
    // Add adds a new particle with the specified parent particle and also passes by reference
    // triangleNormal for adding to m_StructureDirection
    // this is used when a particle attaches to the seed mesh
    void Add(const Vector& p, Vector& triangleNormal, const int parent = -1) {
        const int id = m_Points.size();
        m_Index.insert(make_pair(p.ToBoost(), id));
        m_Points.push_back(p);
        m_Parents.push_back(parent);
        m_JoinAttempts.push_back(0);
        m_BoundingRadius = max(m_BoundingRadius, p.Length() + m_AttractionDistance);

        /////////////////////////////////////////////////////////////////////////
        // Sets the generation number of the particle attached to a seed face
        m_GenerationParent.push_back(1);
        //cout << endl << "parent passed to add as seed";
    // stores the normal of the triangle in the mesh that the seed particle has attached to
        m_StructureDirection.push_back(triangleNormal);
        //cout << endl << m_StructureDirection[m_StructureDirection.size() - 1].X() << endl;
    /////////////////////////////////////////////////////////////////////////



    // update and display progress bar
        ++progressBar;
        progressBar.display();

        // wrap up the progress bar when last particle is placed
        if (id == numParticles) {
            progressBar.done();
        }
    }


    // Nearest returns the index of the particle nearest the specified point
    int Nearest(const Vector& point) const {
        int result = -1;
        m_Index.query(
            boost::geometry::index::nearest(point.ToBoost(), 1),
            boost::make_function_output_iterator([&result](const auto& value) {
                result = value.second;
                }));
        return result;
    }

    vector<int> NearestSmall_n(const Vector& point) const {
        vector<int> result;
        m_IndexSmall.query(
            boost::geometry::index::nearest(point.ToBoost(), knn_slice),
            boost::make_function_output_iterator([&result](const auto& value) {
                result.push_back(value.second);
                }));
        return result;
    }

    // Nearest returns the index of the particle nearest the specified point
    // for the last iterations of RDLA
    Vector NearestSmall(const Vector& point) const {
        //int result = -1;
        Vector result;
        //BoostPoint resultTemp;
        m_IndexSmall.query(
            boost::geometry::index::nearest(point.ToBoost(), 1),
            boost::make_function_output_iterator([&result](const auto& value) {
                //result = value.first;
                result = Vector(boost::geometry::get<0>(value.first),
                    boost::geometry::get<1>(value.first),
                    boost::geometry::get<2>(value.first));
                }));
        return result;
    }

    // Nearest returns the index of the particle nearest the specified point
// for the last iterations of RDLA
    vector<Vector> NearestSmallArray(const Vector& point) const {
        //int result = -1;
        vector<Vector> result;
        //BoostPoint resultTemp;
        m_IndexSmall.query(
            boost::geometry::index::nearest(point.ToBoost(), knn_slice),
            boost::make_function_output_iterator([&result](const auto& value) {
                //result = value.first;
                result.push_back(Vector(boost::geometry::get<0>(value.first),
                    boost::geometry::get<1>(value.first),
                    boost::geometry::get<2>(value.first)));
                }));
        return result;
    }
    

    // DistanceToNearestFace calculates the shortest distance from a particle (p) and nearest face of seed mesh
    // TODO: if needed, consider putting mesh vertices in spatial index, then doing knn search within radius defined by largest edge found in mesh (precomputed)
    double DistanceToNearestFace(const Vector& p, Vector& triangleNormal) {
        double distance = m_AttractionDistance;

        // go over all the shapes in the seed mesh
        for (size_t s = 0; s < seedModelShapes.size(); s++) {

            size_t index_offset = 0;

            // go over all the faces in the shape
            for (size_t f = 0; f < seedModelShapes[s].mesh.num_face_vertices.size(); f++) {
                unsigned int fv = seedModelShapes[s].mesh.num_face_vertices[f];
                vector<Vector> vertices;

                // collect the vertices
                for (size_t v = 0; v < fv; v++) {
                    unsigned int idx = static_cast<unsigned int>(seedModelShapes[s].mesh.indices[index_offset + v].vertex_index);
                    double vx = seedModelAttrib.vertices[3 * idx + 0];
                    double vy = seedModelAttrib.vertices[3 * idx + 1];
                    double vz = seedModelAttrib.vertices[3 * idx + 2];
                    vertices.push_back(Vector(vx, vy, vz));
                }

                index_offset += fv;

                // calculate distance from point to face
                
                double distanceToTriangle = DistanceToTriangle(vertices[0], vertices[1], vertices[2], p);

                ////////////////////////////////////////////////////////////////////////
                // calculates normal of the triangle particle is attaching to, saves the closest triangle 
                // to particle in parameter triangleNormal passed by reference to modify it outside of function and
                // keep updating distance until the smallest distance between point and triangle is found
                if (distanceToTriangle < distance) {
                    distance = distanceToTriangle;

                    // calculate the normal of the closest face to the point
                    // https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
                    // https://math.stackexchange.com/questions/305642/how-to-find-surface-normal-of-a-triangle

                    Vector triangleSide1;
                    Vector triangleSide2;
                    triangleSide1 = vertices[1].operator-(vertices[0]);
                    triangleSide2 = vertices[2].operator-(vertices[0]);

                    // normalise the triangle normal and store in 
                    triangleNormal = (triangleSide1.GetCrossed(triangleSide2)).Normalized();
                }
                
                //distance = min(distance, DistanceToTriangle(vertices[0], vertices[1], vertices[2], p));
            }
        }

        return distance;
    }


    // DistanceToNearestFaceRaw calculates the shortest distance from a particle (p) and nearest face of seed mesh
// TODO: if needed, consider putting mesh vertices in spatial index, then doing knn search within radius defined by largest edge found in mesh (precomputed)
    double DistanceToNearestFaceRaw(const Vector& p) {

        double distance = 100000.0;

        // go over all the shapes in the seed mesh
        for (size_t s2 = 0; s2 < seedModelShapes.size(); s2++) {

            size_t index_offset2 = 0;

            // go over all the faces in the shape
            for (size_t f2 = 0; f2 < seedModelShapes[s2].mesh.num_face_vertices.size(); f2++) {
                unsigned int fv2 = seedModelShapes[s2].mesh.num_face_vertices[f2];
                vector<Vector> vertices;

                // collect the vertices
                for (size_t v2 = 0; v2 < fv2; v2++) {
                    tinyobj::index_t idx2 = seedModelShapes[s2].mesh.indices[index_offset2 + v2];
                    tinyobj::real_t vx2 = seedModelAttrib.vertices[3 * idx2.vertex_index + 0];
                    tinyobj::real_t vy2 = seedModelAttrib.vertices[3 * idx2.vertex_index + 1];
                    tinyobj::real_t vz2 = seedModelAttrib.vertices[3 * idx2.vertex_index + 2];
                    vertices.push_back(Vector(vx2, vy2, vz2));
                }

                index_offset2 += fv2;

                distance = min(distance, DistanceToTriangle(vertices[0], vertices[1], vertices[2], p));
            }
        }
        return distance;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////
 // InitialSeedModelPoint chooses a random face in a random shape (though in our case we only have one shape)
 // then calculates the centre of mass of the face and places the seed point in this location.
 // This function calls the Add() function to output the first point of the DLA simulation
    void InitialSeedModelPoint() {
        //double distance = m_AttractionDistance;
        //cout << static_cast<int>(seedModelShapes.size());
        // choose a random shape in the mesh (there should only be one shape: "0", but this handles the case of multiple shapes
        unsigned int randomShape = RandomInteger(0, static_cast<int>(seedModelShapes.size()) - 1);
        //cout << endl << "randomShape: " << randomShape << endl;

        // choose a random face in a shape
        unsigned int randomFace = RandomInteger(0, static_cast<int>(seedModelShapes[randomShape].mesh.num_face_vertices.size()));
        //cout << endl << "randomFace: " << randomFace << endl;
        // how many vertices does this random face in this random shape have? Should be 3.
        unsigned int fv = seedModelShapes[randomShape].mesh.num_face_vertices[randomFace];
        vector<Vector> vertices;
        // collect the vertices of the randomly chosen face
        for (size_t v = 0; v < fv; v++) {
            tinyobj::index_t idx = seedModelShapes[randomShape].mesh.indices[randomFace + v];
            tinyobj::real_t vx = seedModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 0];
            tinyobj::real_t vy = seedModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 1];
            tinyobj::real_t vz = seedModelAttrib.vertices[static_cast<std::vector<tinyobj::real_t, std::allocator<tinyobj::real_t>>::size_type>(3) * idx.vertex_index + 2];
            vertices.push_back(Vector(vx, vy, vz));
        }
        // find center of mass of the vertices (centre of the triangle)
        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
        for (size_t v = 0; v < vertices.size(); v++) {
            sumX = sumX + vertices[v].X();
            sumY = sumY + vertices[v].Y();
            sumZ = sumZ + vertices[v].Z();
        }
        double centerX = sumX / vertices.size();
        double centerY = sumY / vertices.size();
        double centerZ = sumZ / vertices.size();


        ////////////////////////////////////////////////////////////////////////
        // calculates normal of the triangle particle is attaching to, since
        // seed initial point also needs a normal of the triangle.
        // https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
        // https://math.stackexchange.com/questions/305642/how-to-find-surface-normal-of-a-triangle

        Vector triangleSide1;
        Vector triangleSide2;
        triangleSide1 = vertices[1].operator-(vertices[0]);
        triangleSide2 = vertices[2].operator-(vertices[0]);

        // normalise the triangle normal and store in 
        Vector triangleNormal = (triangleSide1.GetCrossed(triangleSide2)).Normalized();


        Add(Vector(centerX, centerY, centerZ), triangleNormal, -1);
        std::cout << "coordinates of first point: " << centerX << " " << centerY << " " << centerZ << endl;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////




    // RandomStartingPosition returns a random point to start a new particle
    Vector RandomStartingPosition() const {
        const double d = m_BoundingRadius;
        return RandomInUnitSphere().Normalized() * d;
    }

    // ShouldReset returns true if the particle has gone too far away and
    // should be reset to a new random starting position
    bool ShouldReset(const Vector& p) const {
        return p.Length() > m_BoundingRadius * 2;
    }

    // ShouldJoin returns true if the point should attach to the specified
    // parent particle. This is only called when the point is already within
    // the required attraction distance.
    bool ShouldJoin(const Vector& p, const int parent) {
        m_JoinAttempts[parent]++;
        if (m_JoinAttempts[parent] < m_Stubbornness) {
            return false;
        }
        return Random() <= m_Stickiness;
    }

    // PlaceParticle computes the final placement of the particle.
    Vector PlaceParticle(const Vector& p, const int parent) const {
        return Lerp(m_Points[parent], p, m_ParticleSpacing);
    }

    // MotionVector returns a vector specifying the direction that the
    // particle should move for one iteration. The distance that it will move
    // is determined by the algorithm.
    Vector MotionVector(const Vector& p) const {
        return RandomInUnitSphere();
    }

    bool isEqual(double x, double y)
    {
        const double epsilon = 0.00001;
        return abs(x - y) <= epsilon;
    }



    // OutputPointData creates a new point data file and fills it with most current point data
    void OutputPointData(const int iteration) const {
        ///////////////////////////////////////////////////////////////////////////////////////////////////
               // opens file with following name, using the currentDateTime() function above, which prints the current
               // date and time into the file name
        ofstream file;
        //file.open("data/" + seedModelFilename + "_" + filename + "-" + to_string(iteration) + "pts_" + "B" + to_string_with_precision(DefaultStubbornness) +
        //    "K" + to_string_with_precision(DefaultStickiness) + "CA" + to_string_with_precision(DefaultConeAngleDegree) + "CR" + 
        //    to_string_with_precision(DefaultConeReductionRate) + "RD" + to_string_with_precision(DefaultDistanceDiffusionStarts) +
        //    "DATE" + currentDateTime() + ".csv");
        file.open("data/" + filename + "-" + to_string(iteration) + "pts_" + "B" + to_string_with_precision(m_Stubbornness) +
            "K" + to_string_with_precision(m_Stickiness) + "CA" + to_string_with_precision(m_ConeAngleDegree) + "CR" +
            to_string_with_precision(m_ConeReductionRate) + "DD" + to_string_with_precision(m_DistanceDiffusionStarts) +
            "DATE" + currentDateTime() + ".csv");
        if (file.is_open())
            cout << endl << "file opened";
        {
            cout << "OutputPointData function successfully opened file" << endl;
            /*cout << "data/" + seedModelFilename + "_" + filename + "-" + to_string(iteration) + "pts_" + "B" + to_string_with_precision(DefaultStubbornness) +
                "K" + to_string_with_precision(DefaultStickiness) + "CA" + to_string_with_precision(DefaultConeAngleDegree) + "CR" +
                to_string_with_precision(DefaultConeReductionRate) + "RD" + to_string_with_precision(DefaultDistanceDiffusionStarts) +
                "DATE" + currentDateTime() + ".csv";*/
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        for (unsigned int id = 1; id < m_Points.size(); id++) {
            //file << id << "," << m_Parents[id] << "," << m_Points[id].X() << "," << m_Points[id].Y() << "," << m_Points[id].Z() << endl;
            file << m_Points[id].X() << "," << m_Points[id].Y() << "," << m_Points[id].Z() << endl;
        }

        file.close();
    }


    // OutputPointData creates a new point data file and fills it with most current point data
    void OutputDistance(const int iteration) const {
        ///////////////////////////////////////////////////////////////////////////////////////////////////
               // opens file with following name, using the currentDateTime() function above, which prints the current
               // date and time into the file name
        ofstream file;
        file.open("data/" + filename + "-distance_" + to_string(iteration) + "pts_" + "B" + to_string_with_precision(m_Stubbornness) +
            "K" + to_string_with_precision(m_Stickiness) + "CA" + to_string_with_precision(m_ConeAngleDegree) + "CR" +
            to_string_with_precision(m_ConeReductionRate) + "DD" + to_string_with_precision(m_DistanceDiffusionStarts) +
            "DATE" + currentDateTime() + ".csv");
        if (file.is_open())
            cout << endl << "file opened";
        {
            cout << "OutputDistance function successfully opened file" << endl;
            /*cout << "data/" + seedModelFilename + "_" + filename + "-" + to_string(iteration) + "pts_" + "B" + to_string_with_precision(DefaultStubbornness) +
                "K" + to_string_with_precision(DefaultStickiness) + "CA" + to_string_with_precision(DefaultConeAngleDegree) + "CR" +
                to_string_with_precision(DefaultConeReductionRate) + "RD" + to_string_with_precision(DefaultDistanceDiffusionStarts) +
                "DATE" + currentDateTime() + ".csv";*/
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        for (unsigned int id = 0; id < m_Points.size(); id++) {
            //file << id << "," << m_Parents[id] << "," << m_Points[id].X() << "," << m_Points[id].Y() << "," << m_Points[id].Z() << endl;
            file << distanceVector[id] << endl;
        }

        file.close();
    }



    // AddParticle diffuses one new particle and adds it to the model
    void AddParticle(const int iteration) {
        // compute particle starting location
        //Vector p = RandomStartingPosition();
        label2:
        index_error = 0;

        Vector p = RandomInMesh();

        if (iteration > numParticles * attractionBoost) {
            m_AttractionDistance = m_AttractionDistance * 1.5;
            attractionBoost = attractionBoost + 0.01;
        }
            
            

        label1:

        // timer to avoid RDLA getting stuck and to signal that the RDLA is probably full now.
        /*
        const std::chrono::steady_clock::time_point start_time2 = std::chrono::steady_clock::now();
        int offset = 0;
        int attractionOffset = 0;
        //int compoundCounter = 0;
        */
        // do the random walk
        while (true) {
            // get distance to nearest other particle
            /*
            std::chrono::steady_clock::time_point now2 = std::chrono::steady_clock::now();
            auto time_elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(now2 - start_time2).count();

            int time_elapsed_counter = time_elapsed2 - offset;

            // for every 10 seconds elapsed in time_elapsed2
            if ((time_elapsed_counter / 10000) > 1) {
                cout << (time_elapsed2 / 1000) << " seconds elapsed since last point added" << endl;
                offset = offset + 10000;
            }

            // increase attraction distance as time increases. Starts at 2 seconds. Each 2 second that passes
            // the attraction distance is multiplied by 100% with compound interest
            int time_elapsed_attraction_counter = time_elapsed2 - attractionOffset;

            if (time_elapsed_attraction_counter / 1000 > 3) {
                m_AttractionDistance = m_AttractionDistance * 1.5; //* pow(1.3, compoundCounter) ;
                attractionOffset = attractionOffset + 3000;
                //compoundCounter = compoundCounter + 1;
                cout << "Attraction Distance now: " << m_AttractionDistance << endl;
            }



            // if the time for this iteration of AddParticle exceeds 10 seconds, then we quit the RDLA run completely
            // as we have likely finished the RDLA simulation (all slices are full).
            if ((time_elapsed2 / 1000) > 8) {
                // no need for early stop, if we don't want an early stop
                //if (iteration >= numParticles * earlyStop) {
                OutputPointData(iteration);
                OutputDistance(iteration);
                cout << "final printed" << endl;

                for (int i = 0; i < (sizeof(ptsPerSliceCounter) / sizeof(ptsPerSliceCounter[0])); i++) {
                    cout << ptsPerSliceCounter[i] << endl;
                }
                numParticles = 1;
                return;
            }
            */

            int parent;
            double d;

            // if number of iterations exceed the percentage of numParticles defined in
            // smallIndexPerc, then store all the points in the slices that are yet to
            // be full in ptsPerSliceCounter, in the new index m_IndexSmall. This will replace
            // the usual parent and distance d used for searching over for parents to attach to
            if (iteration > smallIndexPerc * numParticles) {
                if (IndexSmallDone == false) {

                    // poitnsSlice was initialised in the inner loop with a dimension that was larger than needed
                    // as the full number of points needed is unknown for all seed meshes. So we need to shrink it
                    // the first resize is just a precaution (though if it didn't work, this would delete data)
                    // the second resize in the loop is for each inner vector, so that it shrinks to the number of 
                    // points designated for each vin
                    //pointsSlice.resize(31);
                    //for (int i = 0; i < ptsPerSlice.size(); ++i) {
                    //    pointsSlice[i].resize(ptsPerSlice[i]);
                    //}

                    // we define iterators for both the outer vector and inner vector (row and col respectively)
                    vector< vector<Vector> >::iterator row;
                    vector<Vector>::iterator col;
                    // these integers are for keeping track of which element of the outer (j) and inner (k) loop
                    // we are in
                    int j = 0;
                    int k = 0;
                    int last_j = 0;
                    // we use the iterators to loop over, starting with the outer loop, the 31 bins
                    for (row = pointsSlice.begin(); row != pointsSlice.end(); row++, j++) {
                        // if the slice is not filled...
                        //cout << j << endl;
                        if (ptsPerSliceCounter[j] < ptsPerSlice[j]) {
                            // print the number of points already in the slice
                            //cout << "counter: " << ptsPerSliceCounter[j] << ", " << "slicesize: " << ptsPerSlice[j] << endl;
                            // lastCol is for printing the size of the bin later on
                            int lastCol = 0;

                            // if the the first element of the inner loop is not full... (this is because we 
                            // are implementing a method for adding the slice immediately below the current
                            // slice in the loop to m_IndexSmall. This is because the lower slice must be searched over
                            // for the candidate point, which is likely to be in the slice above (our current slice) 
                            // that is yet to be filled, to attach to the slice that has vacancies. However, for
                            // the first slice, this doesn't have a previous slice, so this extra condition is needed)
                            if (row != pointsSlice.begin()) {
                                // if the previous outer loop index is just one less than the current (that is, the previous
                                // slice that was not full is immediately below and adjacent to the current slice)...
                                if (last_j == j - 1) {
                                    //cout << "true" << last_j << endl;
                                    //cout << "true" << j << endl;
                                    // we add only the points from the current slice to m_IndexSmall 


                                    for (col = row->begin(); col != row->end(); col++, k++) {
                                        //cout << pointsSlice[j][0].X() << " " << pointsSlice[j][0].Y() << " " << pointsSlice[j][0].Z() << endl;
                                        if (!isEqual((pointsSlice[j][k]).Length(), 0)) {
                                            m_IndexSmall.insert(make_pair(pointsSlice[j][k].ToBoost(), id2));
                                            //cout << pointsSlice[j][k].X() << " " << pointsSlice[j][k].Y() << " " << pointsSlice[j][k].Z() << endl;
                                            id2 = id2 + 1;
                                            lastCol = distance(row->begin(), col);  // Store the current col iterator
                                        }
                                    }
                                    /*
                                    if (lastCol != distance(row->begin(), row->end())) {
                                        //cout << "Last value of col: " << k << endl;
                                        cout << "linked" << endl;
                                        cout << "counter: " << ptsPerSliceCounter[j] << ", " << "slicesize: " << ptsPerSlice[j] << endl;
                                    }
                                    */
                                }
                                // else if the previous outer loop index is not one less than the current...
                                else {
                                    //cout << "false" << last_j << endl;
                                    //cout << "false" << j << endl;
                                    // add both the points from the slice below the current AND the points from
                                    // the current slice. 
                                    // previousInnerVector is for iterating over the outer index value of the previous slice
                                    vector<Vector>& previousInnerVector = *(row - 1);
                                    //cout << !isEqual((pointsSlice[j][k]).Length(), 0) << endl;
                                    for (col = previousInnerVector.begin(); col != previousInnerVector.end(); col++, k++) {
                                        //cout << pointsSlice[j-1][k].X() << " " << pointsSlice[j-1][k].Y() << " " << pointsSlice[j-1][k].Z() << endl;
                                        if (!isEqual((pointsSlice[j - 1][k]).Length(), 0)) {
                                            m_IndexSmall.insert(make_pair(pointsSlice[j - 1][k].ToBoost(), id2));
                                            //cout << pointsSlice[j - 1][k].X() << " " << pointsSlice[j - 1][k].Y() << " " << pointsSlice[j - 1][k].Z() << endl;
                                            id2 = id2 + 1;
                                            lastCol = distance(row->begin(), col);  // Store the current col iterator
                                        }
                                    }
                                    if (lastCol != distance(row->begin(), row->end())) {
                                        //cout << "isolated not full 1" << endl;
                                        //cout << "counter: " << ptsPerSliceCounter[j-1] << ", " << "slicesize: " << ptsPerSlice[j-1] << endl;
                                    }
                                    // reset the k index 
                                    k = 0;
                                    //cout << isEqual((pointsSlice[j][k]).Length(), 0) << endl;
                                    for (col = row->begin(); col != row->end(); col++, k++) {
                                        //cout << pointsSlice[j][k].X() << " " << pointsSlice[j][k].Y() << " " << pointsSlice[j][k].Z() << endl;
                                        if (!isEqual((pointsSlice[j][k]).Length(), 0)) {
                                            m_IndexSmall.insert(make_pair(pointsSlice[j][k].ToBoost(), id2));
                                            //cout << pointsSlice[j][k].X() << " " << pointsSlice[j][k].Y() << " " << pointsSlice[j][k].Z() << endl;
                                            id2 = id2 + 1;
                                            lastCol = distance(row->begin(), col);  // Store the current col iterator
                                        }
                                    }
                                    if (lastCol != distance(row->begin(), row->end())) {
                                        //cout << "isolated not full 2" << endl;
                                        //cout << "counter: " << ptsPerSliceCounter[j] << ", " << "slicesize: " << ptsPerSlice[j] << endl;
                                    }
                                }
                            }
                            // for the first slice (if the first slice is not full, which is unlikely at
                            // a late stage of the RDLA).
                            else {
                                //cout << isEqual((pointsSlice[j][k]).Length(), 0) << endl;
                                for (col = row->begin(); col != row->end(); col++, k++) {
                                    //cout << pointsSlice[j][k].X() << " " << pointsSlice[j][k].Y() << " " << pointsSlice[j][k].Z() << endl;
                                    if (!isEqual((pointsSlice[j][k]).Length(), 0)) {
                                        m_IndexSmall.insert(make_pair(pointsSlice[j][k].ToBoost(), id2));
                                        //cout << pointsSlice[j][k].X() << " " << pointsSlice[j][k].Y() << " " << pointsSlice[j][k].Z() << endl;
                                        id2 = id2 + 1;
                                        lastCol = distance(row->begin(), col);  // Store the current col iterator
                                    }
                                }
                                // Print the last value of col outside the inner loop
                                if (lastCol != distance(row->begin(), row->end())) {
                                    //cout << "first not full" << endl;
                                    cout << "counter: " << ptsPerSliceCounter[j] << ", " << "slicesize: " << ptsPerSlice[j] << endl;
                                }
                            }

                            // integer for comparing the current outer vector index to the index of the outer vector
                            // the next time (ptsPerSliceCounter[j] < ptsPerSlice[j]) is true.



                            // reset the k index 
                            k = 0;
                            last_j = j;
                        }
                    }
                    IndexSmallDone = true;
                    cout << "smallindex done" << endl;
                    cout << "size of m_IndexSmall: " << id2 << endl;
                }

                const Vector pSmall = NearestSmall(p);
                parent = Nearest(pSmall);
                d = p.Distance(m_Points[parent]);
            }
            else {
                parent = Nearest(p);
                d = p.Distance(m_Points[parent]);
            }
        

            // ONLY allow particle to stick when its inside of the base model mesh
            if (IsInsideMesh(p)) {

                // vector for storing the normal of the triangle nearest to the point p
                Vector triangleNormal;

                // check for particle-particle collisions
                if (d < m_AttractionDistance) {

                    // calculate the position of the point, p, m_AttractionDistance from the parent points,
                    // in the direction from parent to p.
                    Vector pRaw = Lerp(m_Points[parent], p, m_ParticleSpacing);

                    // get distance to the nearest seed face
                    const double dfRaw = DistanceToNearestFaceRaw(pRaw);

                    int thicknessRefIndex = findIndex(dfRaw, thicknessRef, 32);

                    // if we get an index error, something has gone wrong with the point, and we can
                    // reset the position of p by giving it a new random position, which will be 
                    // implemented by going to label2 
                    if (index_error == 1) {
                        goto label2;
                    }

                    //cout << thicknessRefIndex << endl;

                    // first condition: if the stubbornness and sticking parameters allow joining with parent
                    // then this cond performs particle placement and add().
                    // originally, this was inverted so the condition was !Shouldjoin(), combined with a continue
                    // statement and no return here. It was used to push the particle away instead of allowing
                    // sticking.
                    /*
                    if (iteration > DLAModePerc * numParticles) {
                        if (ShouldJoin(p, parent)) {

                            //p = Lerp(m_Points[parent], p, m_AttractionDistance + m_MinMoveDistance);
                            p = PlaceParticle(p, parent);
                            Add(p, parent);
                            distanceVector[iteration] = DistanceToNearestFaceRaw(p);
                            return;
                        }
                        // push the particle away if it should not join according to stickiness and stubbornness
                        p = Lerp(m_Points[parent], p, m_AttractionDistance + m_MinMoveDistance);
                        continue;
                    }
                    */
                    // second condition: if random no. matches with laser lichen probability parameters then
                    // allow joining with parent and so this cond is skipped. Otherwise the point is pushed away a bit

                    //cout << thicknessRefIndex << ": " << ptsPerSliceCounter[thicknessRefIndex] << endl;
                    if (!(ptsPerSliceCounter[thicknessRefIndex] <= ptsPerSlice[thicknessRefIndex])) {
                    // push particle away a bit
                        //cout << "counter: " << ptsPerSliceCounter[thicknessRefIndex] << "max: " << ptsPerSlice[thicknessRefIndex] << "index: " << thicknessRefIndex << endl;
                        
                        // if RDLA has exceeded the fastSearchModePerc percentage of the full number of particles
                        // then the RDLA switches into a different mode, where it searches over the knn_slice number of
                        // nearest neighbouring points to p. The indices of the parents closest are placed in array 
                        // parent_n. A for loop goes over the array and does the lerp to see if projecting the point 
                        // until m_AttractionDistance away will allow the point, now called pRawSlice to be in an allowed
                        // slice (not full). So we compare the slice that the pRawSlice is in to how many are already in slice
                        // (ptsPerSliceCounter[thicknessRefIndex]) and the max (ptsPerSlice[thicknessRefIndex]). 
                        // If the pRawSlice is in an allowed slice, then we turn to p and move it 110% m_AttractionDistance 
                        // distance away from the parent point position, in the direction of parent_n to p. Then we continue
                        // so it goes back to top of while loop and walks randomly (if it's still inside the mesh). If none 
                        // of the knn_slice number of points are valid, then we just do the normal lerp with the originally 
                        // found closest point using m_MinMoveDistance etc.
                        
                        if (iteration > fastSearchModePerc * numParticles) {
                            // do search over m_indexSlices

                            // create parent_n, storing the integers of the nearest n number of points to p
                            // in m_IndexSmall. parentSmall_n then holds the coordinates of each of the 
                            // closest points
                            Vector parentVector = Vector();
                            vector<int> parent_n;
                            vector<Vector> parentSmall_n = NearestSmallArray(p);
                            
                            // for each nearest neighbouring knn_slice number of point in m_IndexSmall, find the index 
                            // of the point in m_Index and store it in parent_n
                            for (int m = 0; m < parentSmall_n.size(); m++) {
                                parentVector = parentSmall_n[m];
                                parent_n.push_back(Nearest(parentVector));
                            }

                            // for each index in parent_n, calculate the position of the point if we 
                            // move the point p a m_ParticleSpacing distance from the parent point in
                            // direction of p. 
                            // then check distance of now adjusted p and check which slice it is in 
                            for (int j = 0; j < parent_n.size(); j++) {
                                //cout << j << endl;
                                Vector pRawSlice = Lerp(m_Points[parent_n[j]], p, m_ParticleSpacing);
                                const double dfRawSlice = DistanceToNearestFaceRaw(pRawSlice);
                                unsigned int thicknessRefIndexSlice = findIndex(dfRawSlice, thicknessRef, 32);

                                // if we get an index error, something has gone wrong with the point, and we can
                                // reset the position of p by giving it a new random position, which will be 
                                // implemented by going to label2 
                                if (index_error == 1) {
                                    goto label2;
                                }

                                // If p is in a slice yet to be fully filled, then set p to be in a random direction
                                // around the parent, and go back to the start of the while loop to check if p is 
                                // in a position to attach (check in mesh, check within attach distance and check if 
                                // the random position hasn't moved it to another slice).
                                if ((ptsPerSliceCounter[thicknessRefIndexSlice] <= ptsPerSlice[thicknessRefIndexSlice])) { 
                                    
                                    Vector randDirection = MotionVector(p).Normalized();
                                    //cout << randDirection.X() << ',' << randDirection.Y() << ',' << randDirection.Z() << endl;
                                    p = Lerp(m_Points[parent_n[j]], m_Points[parent_n[j]] + randDirection, m_ParticleSpacing);
                                    //cout << p.X() << ',' << p.Y() << ',' << p.Z() << endl;
                                    goto label1; 
                                }
                                

                                //cout << thicknessRefIndexSlice << endl;

                                //if (!(ptsPerSliceCounter[thicknessRefIndexSlice] >= ptsPerSlice[thicknessRefIndexSlice])) { 
                                //    Vector p = Lerp(m_Points[j], p, m_AttractionDistance * 1.1);
                                //    continue;
                                //}
                            }
                        }
                        // push original particle p away in direction from parent to p since the slice is full
                        // distance is m_AttractionDistance + m_MinMoveDistance
                        p = Lerp(m_Points[parent], p, m_AttractionDistance + m_MinMoveDistance);
                        continue;
                    }
                    
                    
                    
                    // logistic function used to calculate likelihood of attachment style: 
                    // cone based vs RDLA
                    double logisticFunction = 1.0 / (1 + exp(-(dfRaw - m_DistanceDiffusionStarts) / logisticFuncScale));

                    // conditional for whether to use cone-based particle placement or DLA placement
                    if (Random() < (1 - logisticFunction)) {
                        // if random number is smaller than a double that decreases with the generation of parent
                        // place particle using RandomDirectionInCone()
                        // conditional for matching probability of sticking to the histogram of lichen thickness from
                        // the laser scan lichen 

                            //cout << endl <<"direction: " << 
                        //RandomDirectionInCone(parent);
                        //RandomDirectionInCone(parent, 0);
                        // above doesn't work, try another

                        // recalculate distance to the nearest seed face
                        // get point position after random cone direction set
                        Vector point2 = m_Points[parent] + RandomDirectionInCone(parent).Normalized();
                        Vector pRaw2 = Lerp(m_Points[parent], point2, m_ParticleSpacing);

                        // calculate distance of this point from the mesh faces
                        const double dfRaw2 = DistanceToNearestFaceRaw(pRaw2);
                        //cout << dfRaw2 << endl;
                        int thicknessRefIndex2 = findIndex(dfRaw2, thicknessRef, 32);

                        // if we get an index error, something has gone wrong with the point, and we can
                        // reset the position of p by giving it a new random position, which will be 
                        // implemented by going to label2 
                        if (index_error == 1) {
                            goto label2;
                        }

                        // if the slice the point is in is full, then lerp push away, and restart at the 
                        // while(true) {} above 
                        if ((ptsPerSliceCounter[thicknessRefIndex2] > ptsPerSlice[thicknessRefIndex2])) {

                            p = Lerp(m_Points[parent], pRaw2, m_AttractionDistance + m_MinMoveDistance);
                            goto label1;
                        }
                        //cout << " " << RandomDirectionInCone(parent).Normalized().Y();
                        //cout << " " << RandomDirectionInCone(parent).Normalized().Z();

                        // saves the parent index if a particle was placed near it
                        //parent_vector[parentRowIndexCounter][0] = m_Points[parent].X();
                        //parent_vector[parentRowIndexCounter][1] = m_Points[parent].Y();
                        //parent_vector[parentRowIndexCounter][2] = m_Points[parent].Z();

                        //parentRowIndex[parentRowIndexCounter] = parentRow;
                        //++parentRowIndexCounter;
                        //++parentRow;

                        //++faceRow;
                        //cout << "cone" << endl;
                        // add the cone direction point if the slice is not full 
                        Add(pRaw2, parent);
                        // add the point p to the pointsSlice histogram bin thicknessRefIndex2, which has
                        // free space 
                        pointsSlice[thicknessRefIndex2].push_back((m_Points[iteration]));
                        // increment the ptsPerSliceCounter by 1 in the appropriate bin
                        ptsPerSliceCounter[thicknessRefIndex2] = ptsPerSliceCounter[thicknessRefIndex2] + 1;

                        // if the smallIndex is active, then add the point also to the m_IndexSmall index
                        if (IndexSmallDone){
                            m_IndexSmall.insert(make_pair(pRaw2.ToBoost(), id2));
                            id2 = id2 + 1;
                            //cout << m_IndexSmall.size() << endl;
                        }

                        // add the distance of the point pRaw2 to the distance vector, instead of p as pRaw2 was the
                        // point that was attached.
                        distanceVector[iteration] = DistanceToNearestFaceRaw(pRaw2);
                        
                    }

                    else {
                        // else use DLA particle placement
                        //p = PlaceParticle(p, parent);

                        // saves the parent index if a particle was placed near it
                        //parent_vector[parentRowIndexCounter][0] = m_Points[parent].X();
                        //parent_vector[parentRowIndexCounter][1] = m_Points[parent].Y();
                        //parent_vector[parentRowIndexCounter][2] = m_Points[parent].Z();

                        //parentRowIndex[parentRowIndexCounter] = parentRow;
                        //++parentRowIndexCounter;
                        //++parentRow;

                        //++faceRow;
                        //cout << "RDLA" << endl;
                        
                        // calculate point when placed near to face for attaching (replacing PlaceParticle
                        // for transparency).

                        // not needed since duplicate of above
                        //Vector pRaw3 = Lerp(m_Points[parent], p, m_ParticleSpacing);

                        // calculate distance of this point from the mesh faces
                        // not needed, duplicate:
                        //const double dfRaw3 = DistanceToNearestFaceRaw(pRaw3);


                        // find the index of the slice that the point is now in, after it has been brought
                        // closer to the parent for attachment
                        // not needed, duplicate:
                        //int thicknessRefIndex3 = findIndex(dfRaw3, thicknessRef, 32);

                        // if the slice the point is in is full, then lerp push away, and restart at the 
                        // while(true) {} above 
                        // duplicate again, of if slice full comparison:
                        //if ((ptsPerSliceCounter[thicknessRefIndex3] > ptsPerSlice[thicknessRefIndex3])) {
                        //
                        //    p = Lerp(m_Points[parent], pRaw3, m_AttractionDistance + m_MinMoveDistance);
                        //    goto label1;
                        //}

                        // add the point
                        Add(pRaw, parent);
                        pointsSlice[thicknessRefIndex].push_back((m_Points[iteration]));
                        ptsPerSliceCounter[thicknessRefIndex] = ptsPerSliceCounter[thicknessRefIndex] + 1;

                        if (IndexSmallDone) {
                            m_IndexSmall.insert(make_pair(pRaw.ToBoost(), id2));
                            id2 = id2 + 1;
                            //cout << m_IndexSmall.size() << endl;
                        }

                        // update m_StructureDirection to add the cone direction or triangle normal direction 
                        // of the previous parent to the m_Structure. Here it just carries forward the value
                        // of structure direction, in case there is an invocation of the cone attachment mode 
                        // later on in the structure
                        m_StructureDirection.push_back(m_StructureDirection[parent]);

                        distanceVector[iteration] = DistanceToNearestFaceRaw(pRaw);

                    }
                    
                    return;
                }

                // get distance to the nearest seed face
                const double df = DistanceToNearestFace(p, triangleNormal);

                // if the first slice is not full then continue, otherwise ignore this and move on to 
                // moving the point randomly with MotionVector()
                if (ptsPerSliceCounter[0] < ptsPerSlice[0]) {
                    // check for particle-face collisions
                    if (df < m_AttractionDistance) {

                        // calculate the position of the point, p, m_AttractionDistance from the parent points,
                        // in the direction from parent to p.
                        //Vector pRawFace = Lerp(m_Points[parent], p, m_ParticleSpacing);
                        // get distance to the nearest seed face

                        // no need to calculate distance of p from some parent point, since we're attaching to a 
                        // mesh face
                        const double dfRawFace = DistanceToNearestFaceRaw(p);

                        // check if the slice that matches the distance from particle to face (most likely in the bottom-most slice) 
                        // hasn't already been fulled filled. If it's been filled, then push particle away a bit, otherwise, allow
                        // point to stick to face
                        int thicknessRefIndexFace = findIndex(dfRawFace, thicknessRef, 32);

                        // if we get an index error, something has gone wrong with the point, and we can
                        // reset the position of p by giving it a new random position, which will be 
                        // implemented by going to label2 
                        if (index_error == 1) {
                            goto label2;
                        }

                        if (!(ptsPerSliceCounter[thicknessRefIndexFace] <= ptsPerSlice[thicknessRefIndexFace])) {
                            // push particle away a bit
                            //cout << "counter face: " << ptsPerSliceCounter[thicknessRefIndexFace] << "max: " << ptsPerSlice[thicknessRefIndexFace] << "index: " << thicknessRefIndexFace << endl;
                            // p = Lerp(m_Points[parent], p, m_AttractionDistance + m_MinMoveDistance);
                            // move randomly, using the larger of m_MinMoveDistance and df - m_AttractionDistance 
                            // as the distance and a random direction using MotionVector()
                            const double m = max(m_MinMoveDistance, df - m_AttractionDistance);
                            p += MotionVector(p).Normalized() * m;
                            continue;
                        }

                        // add one to the parentRow, to keep track of particles placed
                        //++parentRow;

                        //faceRowIndex[faceRowIndexCounter] = faceRow;
                        //++faceRowIndexCounter;
                        //++faceRow;

                        // saves the particle if it is within AttractionDistance of a face
                        //face_vector[faceRowIndexCounter][0] = p.X();
                        //face_vector[faceRowIndexCounter][1] = p.Y();
                        //face_vector[faceRowIndexCounter][2] = p.Z();

                        // add point pRawFace, which is the closer to the surface at m_ParticleSpacing distance
                        // to the closest mesh face
                        Add(p, triangleNormal, -1);

                        // with the speed up re-indexing of m_IndexSmall, we need to also update the m_IndexSmall
                        // index by adding the point and index to it.
                        if (IndexSmallDone) {
                            m_IndexSmall.insert(make_pair(p.ToBoost(), id2));
                            id2 = id2 + 1;
                            //cout << m_IndexSmall.size() << endl;
                        }

                        //distanceVector[iteration] = df;
                        //cout << m_Points[iteration - 1].X();
                        //cout << pointsSlice[thicknessRefIndexFace][0];
                           
                        // add the point to the 2D vector histogram pointsSlice, which contains coordinates of each point 
                        // that was attached to the slice.
                        pointsSlice[thicknessRefIndexFace].push_back((m_Points[iteration]));
                        // update slice counter to track how many points in a slices
                        ptsPerSliceCounter[thicknessRefIndexFace] = ptsPerSliceCounter[thicknessRefIndexFace] + 1;

                        // update the distace vector with point from pRawFace instead of p as pRawFace is the actual
                        // position of the attached point.
                        distanceVector[iteration] = DistanceToNearestFaceRaw(p);
                        return;
                    }
                }
            }

            // move randomly
            const double m = max(m_MinMoveDistance, d - 0.003);
            p += MotionVector(p).Normalized() * m;

            // check if particle is too far away, reset if so
            if (ShouldReset(p)) {
                p = RandomStartingPosition();
            }
        }
    }




    // m_Index is the spatial index used to accelerate nearest neighbor queries
    Index m_Index;
    // m_IndexSmall is the spatial index used to accelerate nearest neighbor queries
    // for the last RDLA iterations, to speed up each iteration.
    Index m_IndexSmall;

private:
    // m_ParticleSpacing defines the distance between particles that are
    // joined together
    double m_ParticleSpacing;

    // m_AttractionDistance defines how close together particles must be in
    // order to join together
    double m_AttractionDistance;

    // m_MinMoveDistance defines the minimum distance that a particle will move
    // during its random walk
    double m_MinMoveDistance;

    // m_Stubbornness defines how many interactions must occur before a
    // particle will allow another particle to join to it.
    int m_Stubbornness;

    // m_Stickiness defines the probability that a particle will allow another
    // particle to join to it.
    double m_Stickiness;

    // m_BoundingRadius defines the radius of the bounding sphere that bounds
    // all of the particles
    double m_BoundingRadius;

    // m_Points stores the final particle positions
    vector<Vector> m_Points;

    // m_Parents stores the parent IDs of each clustered particle
    vector<int> m_Parents;

    // m_JoinAttempts tracks how many times other particles have attempted to
    // join with each finalized particle
    vector<int> m_JoinAttempts;




    /////////////////////////////////////////////////////////////////////////////////////////
    // NEW PRIVATE MEMBER VARIABLES ADDED BY GAVIN
    // m_generationParent is a private member increment to keep track of how many particles are 
    // already in the structure (structures are DLA point sets with at least one point on a 
    // seed mesh face, each additional point added to the seed point is "adding to the structure").
    // The first point, the seed point on the seed mesh is given m_generationParent = 1.
    vector<int> m_GenerationParent;

    // stores the normal of the triangle that a particle is attaching to or the most recent
    // direction at which a random walker attached to the structure 
    vector<Vector> m_StructureDirection;

    // stores the rate at which the cone width reduces per generation for random vector in cone
    double m_ConeReductionRate;

    // stores cone angle for largest width of the initial non-branch DLA growth
    double m_ConeAngleDegree;

    // stores the rate at which the probability of point attachment changing from cone-based
    // to diffusion based increases (with each generation)
    // DistanceDiffusionStarts is a double between 0.0 and 1.0, where the smaller the number,
    // the less likely the DLA will start earlier on in the structure, and the larger, the
    // more likely DLA will start early in the structure.
    double m_DistanceDiffusionStarts;

};


// create the model in global scope so that parseArgs can configure it
Model model;


// Parses CLI arguments and configures simulation with what is passed
void ParseArgs(int argc, char* argv[]) {
    try {
        cxxopts::Options options(argv[0]);

        options
            .allow_unrecognised_options()
            .add_options()
            ("p,particles", "Number of walker particles", cxxopts::value<int>())
            ("i,input", "Full 3D model filename (.obj)", cxxopts::value<string>())
            ("f,faces", "Seed faces 3D model filename (.obj)", cxxopts::value<string>())
            ("o,output", "Point data output filename", cxxopts::value<string>())
            ("n,interval", "Point data capture interval", cxxopts::value<int>())
            ("s,spacing", "Particle spacing", cxxopts::value<double>())
            ("a,attraction", "Attraction distance", cxxopts::value<double>())
            ("m,move", "Minimum move distance", cxxopts::value<double>())
            ("b,stubbornness", "Stubbornness", cxxopts::value<int>())
            ("k,stickiness", "Stickiness", cxxopts::value<double>())
            ("r,radius", "Initial bounding radius", cxxopts::value<double>())
            ("h,height", "Height of RDLA", cxxopts::value<double>())
            ("c,reduction", "Cone reduction rate", cxxopts::value<double>())
            ("d,angle", "Cone angle degree", cxxopts::value<double>())
            ("l,diffusion", "Distance diffusion starts", cxxopts::value<double>())
            ("w,area", "Area of seed mesh", cxxopts::value<double>())
            ;

        auto result = options.parse(argc, argv);

        if (result.count("particles")) {
            numParticles = result["particles"].as<int>();
            interval = numParticles;
            progressBar.setTotal(numParticles);
        }

        if (result.count("input")) {
            fullModelFilename = result["input"].as<string>();

            // throw error when filename doesn't end in `.obj`
            if (fullModelFilename.substr(fullModelFilename.length() - 4, fullModelFilename.length() - 1).compare(".obj") != 0) {
                cerr << "Base model file must be an OBJ file" << endl;
                exit(1);
            }
        }

        if (result.count("faces")) {
            seedModelFilename = result["faces"].as<string>();

            // throw error when filename doesn't end in `.obj`
            if (seedModelFilename.substr(seedModelFilename.length() - 4, seedModelFilename.length() - 1).compare(".obj") != 0) {
                cerr << "Seed model file must be an OBJ file" << endl;
                exit(1);
            }
        }

        if (result.count("output")) {
            filename = result["output"].as<string>();
        }

        if (result.count("interval")) {
            interval = result["interval"].as<int>();
        }

        if (result.count("spacing")) {
            model.SetParticleSpacing(result["spacing"].as<double>());
        }

        if (result.count("attraction")) {
            model.SetAttractionDistance(result["attraction"].as<double>());
        }

        if (result.count("move")) {
            model.SetMinMoveDistance(result["move"].as<double>());
        }

        if (result.count("stubbornness")) {
            model.SetStubbornness(result["stubbornness"].as<int>());
        }

        if (result.count("stickiness")) {
            model.SetStickiness(result["stickiness"].as<double>());
        }

        if (result.count("radius")) {
            model.SetBoundingRadius(result["radius"].as<double>());
        }
        // to give an option to choose the height in mm of the RDLA growth
        if (result.count("height")) {
            RdlaHeight = result["height"].as<double>();
        }
        // to give an option to choose the cone reduction rate
        if (result.count("reduction")) {
            model.SetConeReductionRate(result["reduction"].as<double>());
        }
        // to give an option to choose the cone angle degree
        if (result.count("angle")) {
            model.SetConeAngleDegree(result["angle"].as<double>());
        }
        // to give an option to choose the height in mm of the RDLA growth
        if (result.count("diffusion")) {
            model.SetDistanceDiffusionStarts(result["diffusion"].as<double>());
        }
        // replacement in case the seed area is calculated incorrectly by the code
        if (result.count("area")) {
            SeedArea = (result["area"].as<double>());
        }
    }
    catch (const cxxopts::OptionException& e) {
        cout << "Error parsing options: " << e.what() << endl;
        exit(1);
    }
}

// LoadFullModel loads the 3D model passed with the `-i` option
void LoadFullModel() {
    string warn;
    string err;

    bool ret = tinyobj::LoadObj(&fullModelAttrib, &fullModelShapes, &fullModelMaterials, &warn, &err, fullModelFilename.c_str());

    if (!warn.empty()) {
        cout << warn << endl;
    }

    if (!err.empty()) {
        cerr << err << endl;
    }

    if (!ret) {
        exit(1);
    }
}

// LoadSeedModel loads the 3D model passed with the `-f` option
void LoadSeedModel() {
    string warn;
    string err;

    bool ret = tinyobj::LoadObj(&seedModelAttrib, &seedModelShapes, &seedModelMaterials, &warn, &err, seedModelFilename.c_str());

    if (!warn.empty()) {
        cout << warn << endl;
    }

    if (!err.empty()) {
        cerr << err << endl;
    }

    if (!ret) {
        exit(1);
    }
}

int main(int argc, char* argv[]) {

    // subtract constant from each element in thicknessRef to account for sheathing performed in Python later
    for (int i = 1; i < 32; ++i) {
        thicknessRef[i] -= sheathThickness;
    }

    // parse the CLI arguments
    ParseArgs(argc, argv);

    // load the full 3D model for inhibiting growth
    if (!fullModelFilename.empty()) {
        LoadFullModel();
    }

    // load the seed faces 3D model
    if (!seedModelFilename.empty()) {
        LoadSeedModel();
    }

    // add seed point at random seed model face (a single point is necessary for now)
    // instead of at origin because that messes up final DLA
    if (!seedModelFilename.empty()) {
        model.InitialSeedModelPoint();
    }
    else
    {
        // add seed point at origin (a single point is necessary for now)
        model.Add(Vector());
    }


    // {
    //     const int n = 3600;
    //     const double r = 1000;
    //     for (int i = 0; i < n; i++) {
    //         const double t = (double)i / n;
    //         const double a = t * 2 * M_PI;
    //         const double x = cos(a) * r;
    //         const double y = sin(a) * r;
    //         model.Add(Vector(x, y, 0));
    //     }
    // }




    // redefine numParticles so that it converts between RDLA height(user input as - h height e.g. - h 0.001)
    // the numParticles will then adjust to the varying area 
    // temporarily removed to calculate the number of points against height 
    // NEED TO CHANGE THE AREA DIVISION, WHICH IS 0.0348274, DIVIDE IT BY THAT TO GET SCALED NUMBER OF POINTS CONVERSION
    // CHANGE TO NUMPARTICLES CALCULATION:
    // we will now do area of seed divided by ((0.01^2)*pi)*200. This is the number of points (200) in the z-collapsed 
    // cylinder of points in circle of radius 1 cm. 
    //numParticles = int((CalculateTotalFaceArea() / (pow(0.01, 2) * pi)) * 200);
    numParticles = int((SeedArea / (pow(0.01, 2) * pi)) * 200);
    //numParticles = int(((RdlaHeight - (2.36 * pow(10, -3))) / (3.64 * pow(10, -7))) * (CalculateTotalFaceArea() / 0.0461113)); // 0.0461113 is area of the carving F596, which the RDLA height-to-points conversion is based on 


    //cout << totalProbability << endl;
    // loop for assigning the maximum number of points allowed per slice, as defined by the histogram created from the
    // lichen laser scan. The max number is calculated by normalising all probabilities (which was already converted from
    // counts of histogram, so this is an unnecessary extra step, could just normalise counts). Then, we divide each probability
    // with the total probability and then scale to the max number of particles, based on the height of lichen wanted. 
    int numParticlesTemp = 0;
    
    for (int i = 0; i < size(thicknessProbability); i++) {
        ptsPerSlice[i] = floor((thicknessProbability[i] / totalProbability) * numParticles); // (meshArea * thicknessRef[i])/pow(m_ParticleSpacing,3)
        numParticlesTemp = numParticlesTemp + ptsPerSlice[i];
        cout << ptsPerSlice[i] << endl;
    }

    numParticles = numParticlesTemp;

    cout << numParticles << endl;

    

    // run diffusion-limited aggregation
    //cout << CalculateTotalFaceArea() << endl;

    distanceVector.resize(numParticles);

    // define the dimensions of the pointSlice, storing the coordinates of the attached points
    //int maxElementValue = *std::max_element(ptsPerSlice.begin(), ptsPerSlice.end());

    //vector<vector<Vector>> pointsSlice;

    cout << "Total area of mesh: " << SeedArea << endl;

    //CalculateTotalFaceAreaPrint();



    for (int i = 1; i <= numParticles; i++) {
        model.AddParticle(i);

        if (i == numParticles) {
            // no need for early stop, if we don't want an early stop
            //if (iteration >= numParticles * earlyStop) {
            model.OutputPointData(i);
            model.OutputDistance(i);
            cout << "final printed" << endl;

            for (int j = 0; j < (sizeof(ptsPerSliceCounter) / sizeof(ptsPerSliceCounter[0])); j++) {
                cout << ptsPerSliceCounter[j] << endl;
            }
            //numParticles = 1;
            //return;
        }
        // output current point data based on interval
        /*
        if (i % 1000 == 0) {
            model.OutputPointData(i);
            cout << "printed";
            //lastOutputIteration = i;

            //cout << endl << "Gamma value (for area and height to point number conversion): ";
            //GammaCalculate();

            //cout << endl << "Total area of mesh: " << CalculateTotalFaceArea();

            // for printing out the parent points to a file to check 
            // why there is an issue with particle being placed too far 
            // from its parent 

            //file.open("data/" + filename + "_parent_" + to_string(i) + "pts_" + currentDateTime() + ".csv", ios::app);

            //for (size_t k = 0; k < parentRowIndexCounter; ++k)
            //{
            //    //cout << endl << k << endl;
            //    file << parentRowIndex[k] << ",";
            //    for (size_t j = 0; j < 3; ++j) {
            //        if (j < (3 - 1)) {
            //            file << parent_vector[k][j] << ",";
            //        }
            //        else if (j == (3 - 1)) {
            //            file << parent_vector[k][j] << "\n";
            //        }
            //    }
            //}
            //cout << "parent file successfully opened" << endl;
            //file.close();


            //file.open("data/" + filename + "_face_" + to_string(i) + "pts_" + currentDateTime() + ".csv", ios::app);

            //for (size_t k = 0; k < faceRowIndexCounter; ++k)
            //{
            //    //cout << endl << k << endl;
            //    file << faceRowIndex[k] << ",";
            //    for (size_t j = 0; j < 3; ++j) {
            //        if (j < (3 - 1)) {
            //            file << face_vector[k][j] << ",";
            //        }
            //        else if (j == (3 - 1)) {
            //            file << face_vector[k][j] << "\n";
            //        }
            //    }
            //}
            //cout << "face file successfully opened" << endl;
            //file.close();

        }*/
        

    }

    return 0;
}
