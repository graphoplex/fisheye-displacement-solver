function median_discrepancy = mini(matchedPoints1,d_observed, inputs)
%MINI function to minimize via fmincon
%with inputs = dx, dy, dtheta, x0,y0,f, h, k1,k2,k3,k4
%with X = x, y, dtheta, x0,y0,f, h, k1,k2,k3,k4
r = @(X) X(1)^2/X(7)^2 + X(2)^2/X(7)^2;
alpha = @(X) atan(r(X(1:7)));
alphad = @(X) alpha(X(1:7))*(1+X(8)*alpha(X(1:7))^2+X(9)*X(1)^4+X(10)*alpha(X(1:7))^6+X(11)*alpha(X(1:7))^8);
xp = @(X)(alphad(X)/r(X))*(X(1)/X(7));
yp = @(X)(alphad(X)/r(X))*(X(2)/X(7));
u = @(X) X(6)*(X(4)-yp(X)*cos(X(3))-xp(X)*sin(X(3)));
v = @(X) X(6)*(X(5)-xp(X)*cos(X(3))+yp(X)*sin(X(3)));

%reverse
%with X = u,v,x0,y0,f
r_yp = @(X) X(3)-X(1)/X(5);
r_xp = @(X) X(4)-X(2)/X(5);

%simulated displacement
%with X = x, y, dtheta, x0,y0,f, h, k1,k2,k3,k4,u,v
newxp = @(X)xp(X(1:11)) + r_xp([X(12),X(13),X(4:6)]);
newyp = @(X)yp(X(1:11)) + r_yp([X(12),X(13),X(4:6)]);
newu = @(X) X(6)*(X(4)-yp(X)*cos(X(3))-xp(X)*sin(X(3)));
newv = @(X) X(6)*(X(5)-xp(X)*cos(X(3))+yp(X)*sin(X(3)));

inputs
discrepancies = zeros(size(matchedPoints1));
for i = 1:size(matchedPoints1)
    feature = matchedPoints1(i);
    ufeature = feature.Location(1);
    vfeature = feature.Location(2);
    rx = r_xp([ufeature, vfeature, inputs(4:6)]);
    ry = r_yp([ufeature, vfeature, inputs(4:6)]);
    syms x  y;
    assume(10*x,'integer');
    assume(10*y,'integer');
    sols = vpasolve(rx == x* inputs(7)*alphad([x,y,inputs(3:11)])/(x^2+y^2),ry == y* inputs(7)*alphad([x,y,inputs(3:11)])/(x^2+y^2),x, y);
    newx = sols.x+inputs(1);
    newy = sols.y+inputs(2);
    newufeature = newu([newx,newy,inputs(3:11)]);
    newvfeature = newv([newx,newy,inputs(3:11)]);
    du_simulated = newufeature - ufeature;
    dv_simulated = newvfeature - vfeature;
    discrepancies(i) = sqrt((du_simulated-d_observed(i,1))^2+(dv_simulated-d_observed(i,2))^2);
end
median_discrepancy = sum(discrepancies)      
end

