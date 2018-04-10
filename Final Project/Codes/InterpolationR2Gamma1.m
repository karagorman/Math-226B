% Math 226B - Final Project
% Part 4a
% Function to compute interpolation matrix I_R2Gamma1

function v_Gamma1 = InterpolationR2Gamma1(m,EnumR1,EnumR2,v_R2,tot_R1_pts,...
    tot_Sigma1_pts,totalR1pts)

h = 1/(2^m);
v_Gamma1 = zeros(totalR1pts - (tot_R1_pts + tot_Sigma1_pts),1);
k = 1;
for i = tot_R1_pts + tot_Sigma1_pts + 1:totalR1pts
    [x1,y1] = find(EnumR1 == i);
    x2 = x1 + 1;
    y2 = y1 + 1;
    if (EnumR2(x1,y1+2) == 0)
        y1 = y1 - 1;
        y2 = y1 - 1;
    elseif (EnumR2(x1-1,y1) == 0)
        y1 = y1 + 1;
        y2 = y1 + 1;
    else
        y2 = y1 + 1;
    end
    RHS = [v_R2(EnumR2(x1,y1)); v_R2(EnumR2(x1,y2)); v_R2(EnumR2(x2,y1));...
        v_R2(EnumR2(x2,y2))];
    x1 = x1*h - 2*h; x2 = x2*h - 2*h; y1 = y1*h - 2*h; y2 = y2*h - 2*h; 
    Int_Mat1 = [1,x1,y1,x1*y1; 1,x1,y2,x1*y2; 1,x2,y1,x2*y1; 1,x2,y2,x2*y2];
    c = Int_Mat1\RHS;
    v_Gamma1(k) = c(1) + c(2)*x1 + c(3)*y1 + c(4)*x1*y1;
    k = k+1;
end