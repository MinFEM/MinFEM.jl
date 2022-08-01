//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0.0};
//+
Physical Surface("domain", 10001) = {1};
//+
Physical Line("top", 1001) = {3};
//+
Physical Line("right", 1002) = {2};
//+
Physical Line("bottom", 1003) = {1};
//+
Physical Line("left", 1004) = {4};
//+
