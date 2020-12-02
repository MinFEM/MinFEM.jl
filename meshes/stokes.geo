//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-7, -3, 0, 14, 6, 0};
//+
Disk(3) = {0, 0, 0, 0.5, 0.5};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Delete; }

// Interior volume
Physical Surface(10001) = {1};

// Outer boundary
Physical Line(1001) = {6,9};

// Inlet boundary
Physical Line(1003) = {7};

// Obstacle boundary
Physical Line(1002) = {5};

// Outlet boundary
Physical Line(1004) = {8};

// Nose
Physical Point(101) = {1};

Field[1] = Attractor;
Field[1].EdgesList = {5};
Field[1].NNodesByEdge = 2000;
Field[3] = MathEval;
Field[3].F = "min(log(F1+1.1), 0.5)";
Background Field = 3;
