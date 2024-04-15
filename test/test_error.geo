// Square mesh that does not contain a domain
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Physical Curve("bottom", 1001) = {1};
//+
Physical Curve("top", 1002) = {3};
//+
Physical Curve("left", 1003) = {4};
//+
Physical Curve("right", 1004) = {2};
