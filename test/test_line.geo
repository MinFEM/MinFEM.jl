// Gmsh project created on Fri Jan 07 15:56:13 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, -0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("left", 1001) = {1};
//+
Physical Point("right", 1002) = {2};
//+
Physical Curve("domain", 10001) = {1};
