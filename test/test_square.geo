// Gmsh project created on Wed Jan 12 15:59:38 2022
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Curve("bottom", 1001) = {1};
//+
Physical Curve("top", 1002) = {3};
//+
Physical Curve("left", 1003) = {4};
//+
Physical Curve("right", 1004) = {2};
//+
Physical Surface("domain", 10001) = {1};
