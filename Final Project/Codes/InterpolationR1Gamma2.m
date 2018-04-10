% Math 226B - Final Project
% Part 4a
% Function to compute interpolation matrix I_R1Gamma2

function v_Gamma2 = InterpolationR1Gamma2(m,EnumR1,EnumR2,v_R1,tot_R2_pts,...
    tot_Sigma2_pts,totalR2pts)

h = 1/(2^m);
v_Gamma2 = zeros(totalR2pts - (tot_R2_pts + tot_Sigma2_pts),1);
k = 1;
for i = tot_R2_pts + tot_Sigma2_pts + 1:totalR2pts
    [x1,y1] = find(EnumR2 == i);
    y2 = y1 - 1;
    if (EnumR1(x1,y1-2) == 0);
        x1 = x1 - 1;
        x2 = x1 - 1;
        
        if (EnumR1(x1,y2) > length(v_R1))
            y2 = y1 + 1;
        end
          
    elseif (EnumR1(x1+1,y1) == 0)
        x1 = x1 - 1;
        x2 = x1 - 1;
        y2 = y1 + 1;
    else
        x2 = x1 - 1;
    end
    
    RHS = [v_R1(EnumR1(x1,y1)); v_R1(EnumR1(x1,y2)); v_R1(EnumR1(x2,y1));...
        v_R1(EnumR1(x2,y2))];
    x1 = x1*h - 2*h; x2 = x2*h - 2*h; y1 = y1*h - 2*h; y2 = y2*h - 2*h;
    Int_Mat2 = [1,x1,y1,x1*y1; 1,x1,y2,x1*y2; 1,x2,y1,x2*y1; 1,x2,y2,x2*y2];
    c = Int_Mat2\RHS;
    v_Gamma2(k) = c(1) + c(2)*x1 + c(3)*y1 + c(4)*x1*y1;;
    k = k + 1;
end