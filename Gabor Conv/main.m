clc
clear
Ib=imread('elephant.png');
Ibd = im2double(Ib);
%% Q1
g1 = gabor2(3,0.1,90,0.75,90);
Q1 = conv2(Ibd, g1, 'valid');
q11 = Q1(99, 281);
q12 = Q1(221, 54);
%% Q2
g2 = gabor2(3,0.1,90,0.75,90);
g3 = gabor2(3,0.1,90,0.75,0);
Q21 = conv2(Ibd, g2, 'valid');
Q22 = conv2(Ibd, g3, 'valid');
Q2 = sqrt(Q21.^2 + Q22.^2);
q21 = Q2(253, 165);
q22 = Q2(169, 37);
%% Q3
Q3 = 0;
for orientation = 0:15:179
    g4 = gabor2(3,0.1,orientation,0.75,90);
    g5 = gabor2(3,0.1,orientation,0.75,0);
    Q31 = conv2(Ibd, g4, 'valid');
    Q32 = conv2(Ibd, g5, 'valid');
    Q = sqrt((Q31.^2) + (Q32.^2));
    Q3 = max(Q, Q3);
end
q31 = Q3(472, 315);
q32 = Q3(284, 252);