//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {2, 1, 0, 1.0};
//+
Point(4) = {3, 1, 0, 1.0};
//+
Point(5) = {3, 2, 0, 1.0};
//+
Point(6) = {1, 2, 0, 1.0};
//+
Point(7) = {1, 1, 0, 1.0};
//+
Point(8) = {0, 1, 0, 1.0};
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
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Curve Loop(1) = {6, 7, 8, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom_left", 1001) = {1};
//+
Physical Curve("right_lower", 1002) = {2};
//+
Physical Curve("bottom_right", 1003) = {3};
//+
Physical Curve("right_upper", 1004) = {4};
//+
Physical Curve("top_right", 1005) = {5};
//+
Physical Curve("left_upper", 1006) = {6};
//+
Physical Curve("top_left", 1007) = {7};
//+
Physical Curve("left_lower", 1008) = {8};
//+
Physical Surface("domain", 10001) = {1};
