%Parameters to optimize


% Feature detection and matching
im1=rgb2gray(imread('o1.jpg'));
im2=rgb2gray(imread('o2.jpg'));

points1 = detectHarrisFeatures(im1,'MinQuality',0.003);
points2 = detectHarrisFeatures(im2,'MinQuality',0.003);

[features1,valid_points1] = extractFeatures(im1,points1);
[features2,valid_points2] = extractFeatures(im2,points2);

indexPairs = matchFeatures(features1,features2);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
displacements_observed = matchedPoints2.Location - matchedPoints1.Location;


figure; showMatchedFeatures(im1,im2,matchedPoints1,matchedPoints2);

%___Optimization w/ constraints____
%dx, dy, dtheta, x0,y0,f, h, k1,k2,k3,k4
x0 = [0.,0.,0.,0.,0.,0.001,1.,0.,0.,0.,0.];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-1,-1,-0.1,0 ,0 ,0  ,0,0  ,0  ,0  ,0];
ub = [1 ,1 ,0.1 ,0 ,0 ,0.5,5,100,100,100,100];

f = @(X) mini(matchedPoints1, displacements_observed, X);
displacement = fmincon(f, x0,A,b,Aeq,beq,lb,ub)

