//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0.7, 0, 1, 0.3, 0};
//+
Rectangle(2) = {0, 0.1, 0, 0.8, 0.3, 0};
//+
Rectangle(3) = {-0, -0.5, 0, 0.3, 1.5, 0};
//+
Rectangle(4) = {1.2, -0.5, 0, 0.3, 1.5, 0};
//+
Rectangle(5) = {1.2, 0.1, -0, 0.8, 0.3, 0};
//+
Rectangle(6) = {1.2, 0.7, -0, 1, 0.3, 0};
//+
Rectangle(7) = {1.2, -0.5, -0, 1, 0.3, 0};
//+
Rectangle(8) = {2.4, -0.5, -0, 0.3, 1.5, 0};
//+
Rectangle(9) = {3.7, -0.5, -0, 0.3, 1.5, 0};
//+
Point(37) = {3.2, -0.1, 0, 1.0};
//+
Point(38) = {3.2, -0.5, 0, 1.0};
//+
Point(39) = {2.7, 0.6, 0, 1.0};
//+
Line(37) = {31, 37};
//+
Line(38) = {37, 38};
//+
Line(39) = {38, 39};
//+
Point(40) = {3.7, 0.6, 0, 1.0};
//+
Line(40) = {36, 37};
//+
Line(41) = {38, 40};
//+
BooleanUnion{ Surface{3}; Delete; }{ Surface{2}; Surface{1}; Delete; }
//+
BooleanUnion{ Surface{7}; Delete; }{ Surface{5}; Surface{4}; Surface{6}; Delete; }
//+
Line(86) = {39, 31};
//+
Line(87) = {36, 40};
//+
Curve Loop(24) = {39, 86, 37, 38};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {41, -87, 40, 38};
//+
Curve Loop(26) = {36, 33, 34, 35};
//+
Curve Loop(27) = {41, -87, 40, 38};
//+
Plane Surface(25) = {27};
//+
BooleanUnion{ Surface{8}; Delete; }{ Surface{24}; Surface{25}; Surface{9}; Delete; }
//+
Translate {-6.3, 0, 0} {
  Duplicata { Surface{8}; Surface{24}; Surface{25}; Surface{9}; }
}
//+
Point(115) = {-0.3, -0.5, 0, 1.0};
//+
Point(116) = {-0.6, -0.5, 0, 1.0};
//+
Point(117) = {-1, 0, 0, 1.0};
//+
Point(118) = {-1.3, -0.3, 0, 1.0};
//+
Point(119) = {-1.3, 0.0, -0, 1.0};
//+
Point(120) = {-1, 0.3, -0, 1.0};
//+
Point(121) = {-0.5, 0.2, 0, 1.0};
//+
Point(122) = {-0.3, -0.2, 0, 1.0};
//+
Point(123) = {-0.6, -0.1, 0, 1.0};
//+
BSpline(132) = {116, 116, 123, 117, 118, 118, 119, 119, 120, 121, 122, 115, 115, 115, 116};
//+
Rectangle(30) = {-1.4, -0.5, 0, 0.3, 0.8, 0};
//+
Rectangle(31) = {-2, -0.5, 0, 0.3, 0.8, 0};
//+
Rectangle(32) = {-2, 0.4, 0, 0.3, 0.3, 0};
//+
Curve Loop(39) = {132};
//+
Plane Surface(33) = {39};
//+
BooleanFragments{ Surface{33}; Delete; }{ Surface{30}; Delete; }
//+
Rectangle(36) = {-4.1, -0.75, 0, 8.4, 2, 0};
//+
BooleanDifference{ Surface{36}; Delete; }{ Surface{26}; Surface{27}; Surface{28}; Surface{29}; Surface{31}; Surface{32}; Surface{35}; Surface{33}; Surface{34}; Surface{10}; Surface{11}; Surface{12}; Surface{13}; Surface{15}; Surface{14}; Surface{23}; Surface{22}; Surface{21}; Surface{18}; Surface{19}; Surface{20}; Surface{16}; Surface{17}; Surface{8}; Surface{24}; Surface{25}; Surface{9}; Delete; }
//+
Physical Curve("Outer Boundary", 1001) = {155, 157, 156, 154};
//+
Physical Curve("Inner Boundary", 1002) = {118, 164, 165, 158, 159, 114, 160, 163, 162, 128, 127, 161, 140, 139, 141, 144, 143, 142, 138, 137, 149, 152, 147, 145, 150, 151, 153, 42, 46, 49, 52, 53, 58, 59, 60, 51, 55, 56, 57, 44, 45, 67, 64, 61, 75, 68, 77, 80, 81, 83, 84, 85, 79, 72, 73, 74, 76, 65, 66, 91, 90, 37, 39, 89, 88, 41, 40, 94, 93, 95, 92};
//+
Physical Surface("Volume", 10001) = {36};
//+
Characteristic Length {136} = 0.3;
//+
/*Field[1] = Attractor;*/
/*//+*/
/*Field[1].EdgesList = {132};*/
/*//+*/
/*Field[2] = Threshold;*/
/*//+*/
/*Field[1].NNodesByEdge = 100;*/
/*//+*/
/*Field[2].DistMax = 0.1;*/
/*//+*/
/*Field[2].DistMin = 0;*/
/*//+*/
/*Background Field = 2;*/
/*//+*/
/*Field[2].IField = 1;*/
