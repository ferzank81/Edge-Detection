function [heat_av ] = heattrans(maskGray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved.
% This work should only be used for nonprofit purposes.
% Please cite the paper when you use this code:
% Katýrcýoðlu, F. (2020). Edge detection method based on heat conduction matrix for infrared images.
% Optical Engineering, 59(9), 093103.
%
% AUTHOR:
%     Ferzan Katýrcýoðlu,Duzce University, TURKEY.
%     email:katirciogluferzan@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=9;
kat=0;

e=0.0001;
mask(1,1)=maskGray(1);
mask(1,2)=maskGray(2);
mask(1,3)=maskGray(3);
mask(2,1)=maskGray(4);
mask(2,2)=maskGray(9);
mask(2,3)=maskGray(5);
mask(3,1)=maskGray(6);
mask(3,2)=maskGray(7);
mask(3,3)=maskGray(8);



[max_ust1 sira11]=max(mask(1,:));
[max_ust2 sira12]=max(mask(3,:));
if(max_ust1>max_ust2);
    max_ust=max_ust1;
    sira1=sira11;
else
    max_ust=max_ust2;
    sira1=sira12;
end
[min_alt1 sira21]=min(mask(1,:));
[min_alt2 sira22]=min(mask(3,:));
if(min_alt1<min_alt2);
    min_alt=min_alt1;
    sira2=sira21;
else
    min_alt=min_alt2;
    sira2=sira22;
end

% kat1=(max_ust-maskGray(9));
% kat2=(min_alt-maskGray(9));
% kat=(kat1+kat2)/(64);
% kat=(max_ust+min_alt+maskGray(9))/(255);
kat=(maskGray(9)-mean(maskGray(:)))/4;

heat_av=0;



%  Determining the L value in heat transfer
if (sira1==2)
    L1=1;
    A1=14;
else
    L1=1.4142;
    A1=9;
end
if (sira2==2 )
    L2=1;
    A2=14;
else
    L2=1.4142;
    A2=9;
end
L=L1+L2;
A=(A1+A2)/2;

% Calculation of heat in the mask according to the center pixel

 heat_av=kat*A*[(max_ust-min_alt)]/ L;
 

end

