function [t]=RotateMatrixToEuler(RotateA)
%RotateA is Rotate Matrix. 
%output: t has 3 elements, [x,y,z] rotation
%Bo Li 12-20-2020
theta=asind(RotateA(3,1));
if abs(theta-90)>0.01
    if (RotateA(3,3)*cosd(theta)>0)
        fai=atand(-RotateA(3,2)/RotateA(3,3));
    else
        fai=atand(-RotateA(3,2)/RotateA(3,3))+180;
    end
    if (RotateA(1,1)*cosd(theta)>0)
        posai=atand(-RotateA(2,1)/RotateA(1,1));
    else
        posai=atand(-RotateA(2,1)/RotateA(3,3))+180;
    end
        t=[fai, theta, posai];
else
    t=[-1, theta, -1];
end
    
end