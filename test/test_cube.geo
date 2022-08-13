// Gmsh project created on Fri Jan 07 11:24:20 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {0, 0, 1, 1.0};
//+
Point(6) = {0, 1, 1, 1.0};
//+
Point(7) = {1, 1, 1, 1.0};
//+
Point(8) = {1, 0, 1, 1.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 8};
//+
Line(3) = {8, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {6, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 7};
//+
Line(8) = {7, 6};
//+
Line(9) = {5, 6};
//+
Line(10) = {8, 7};
//+
Line(11) = {1, 2};
//+
Line(12) = {4, 3};
//+
Curve Loop(1) = {9, -8, -10, -2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 5, -11, 1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 6, -12, 4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 7, -10, 3};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, 3, 4, 1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {5, 6, 7, 8};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 2, 6, 3, 4, 5};
//+
Volume(1) = {1};
//+
Physical Surface("front", 1001) = {1};
//+
Physical Surface("right", 1002) = {4};
//+
Physical Surface("back", 1003) = {3};
//+
Physical Surface("left", 1004) = {2};
//+
Physical Surface("bottom", 1005) = {5};
//+
Physical Surface("top", 1006) = {6};
//+
Physical Volume("domain", 10001) = {1};
