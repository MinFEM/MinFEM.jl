//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {3, 0, 0, 1.0};
//+
Point(3) = {3, 2, 0, 1.0};
//+
Point(4) = {2, 2, 0, 1.0};
//+
Point(5) = {2, 1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {5, 6, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 1001) = {1};
//+
Physical Curve("right", 1002) = {2};
//+
Physical Curve("top_right", 1003) = {3};
//+
Physical Curve("left_upper", 1004) = {4};
//+
Physical Curve("top_left", 1005) = {5};
//+
Physical Curve("left_lower", 1006) = {6};
//+
Physical Surface("domain", 10001) = {1};
