
%% 总结
% 1. 无特殊说明， 一般卷积图像参数选择'same'
% 2. Ensure the size of each mask is sufficient to accurately represent 
%    the Gaussian. 意味着产生的高斯mask大小为sigma的6倍（见lecture）
% 3. 保留两位小数， roundn函数可以用来保留小数：roundn(数值， -2)
% 4. 注意用的哪个图像
% 5. same 与 valid 的区别
%% Q1：Smoothing using Box Masks
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('boxes.pgm');

Ia_gray = rgb2gray(Ia); 

Ia_grayd = im2double(Ia_gray); %转换成double类型
Icd = im2double(Ic); %转换成double类型
Ibd = im2double(Ib);

Box5 = fspecial('average', [5 5]);
Box25 = fspecial('average', [25 25]);
convBox5_Ia_grayd = conv2(Ia_grayd, Box5, 'same');
convBox5_Icd = conv2(Icd, Box5, 'same');
convBox25_Ia_grayd = conv2(Ia_grayd, Box25, 'same');
convBox25_Icd = conv2(Icd, Box25, 'same');

figure(1)
subplot(2, 2, 1), imshow(convBox5_Ia_grayd)
subplot(2, 2, 2), imshow(convBox5_Icd)
subplot(2, 2, 3), imshow(convBox25_Ia_grayd)
subplot(2, 2, 4), imshow(convBox25_Icd)
%% Q2: Smoothing using Gaussian Masks
G1_5 = fspecial('gaussian', [9 9], 1.5);
G10 = fspecial('gaussian', [60 60], 10);

convG1_5_Ia_grayd = conv2(Ia_grayd, G1_5, 'same');
convG1_5_Icd = conv2(Icd, G1_5, 'same');
convG10_Ia_grayd = conv2(Ia_grayd, G10, 'same');
convG10_Icd = conv2(Icd, G10, 'same');

figure(2)
subplot(2, 2, 1), imshow(convG1_5_Ia_grayd)
subplot(2, 2, 2), imshow(convG1_5_Icd)
subplot(2, 2, 3), imshow(convG10_Ia_grayd)
subplot(2, 2, 4), imshow(convG10_Icd)
%% Q3: Differencing with First-Derivative Masks
y_fdmask = [-1; 1];
x_fdmask = [-1, 1];

conv_y_fdmask_Ibd = conv2(Ibd, y_fdmask, 'valid');
conv_x_fdmask_Ibd = conv2(Ibd, x_fdmask, 'valid');
%% Q4: Differencing with Laplacian Masks
L1 = -1 * ones(3, 3); L1(2,2) = 8;
L2 = L1 / 8; %标准化后的L mask

conv_L2_Ia_grayd = conv2(Ia_grayd, L2, 'same');
conv_L2_Icd = conv2(Icd, L2, 'same');
Q4 = conv_L2_Icd(22, 41:45);
%% ii
% Sobel mask
Sobel1 = fspecial('Sobel');

conv_Sobel1_Icd = conv2(Icd, Sobel1);
conv_Sobel1T_Icd = conv2(Icd, Sobel1');

figure
subplot(1, 2, 1), imshow(conv_Sobel1_Icd);
subplot(1, 2, 2), imshow(conv_Sobel1T_Icd);

% Prewitt mask
Prewitt1 = fspecial('Prewitt');

conv_Prewitt1_Icd = conv2(Icd, Prewitt1);
conv_Prewitt1T_Icd = conv2(Icd, Prewitt1');

figure
subplot(1, 2, 1), imshow(conv_Prewitt1_Icd);
subplot(1, 2, 2), imshow(conv_Prewitt1T_Icd);
%% Q5: Edge Detection with Gaussian Derivative Masks
G2_5 = fspecial('gaussian', [15, 15], 2.5);
GDmask_x = conv2(G2_5, x_fdmask, 'valid');
GDmask_y = conv2(G2_5, y_fdmask, 'valid');

figure
subplot(2, 2, 1), mesh(GDmask_x);
subplot(2, 2, 2), mesh(GDmask_y);
conv_GDmask_x_Icd = conv2(Icd, GDmask_x, 'same'); %这里要尝试一下对比valid
conv_GDmask_y_Icd = conv2(Icd, GDmask_y, 'same');
subplot(2, 2, 3), imshow(conv_GDmask_x_Icd);
subplot(2, 2, 4), imshow(conv_GDmask_y_Icd);
Q_L2 = sqrt(conv_GDmask_x_Icd.^2+conv_GDmask_y_Icd.^2);
figure, imshow(Q_L2)
Q5 = conv_GDmask_x_Icd(22,41:45);
%% Q6: Edge Detection with Laplacian of Gaussian (LoG) Masks
g1 = fspecial('gaussian', 13, 1.5); %高斯算子足够大，最终的LoG基本一样
g2 = fspecial('gaussian', 45, 5);
g1l = conv2(g1, L2, 'valid'); %这里用valid和same的结果近似，取两位小数后一样
g2l = conv2(g2, L2, 'valid');

figure
subplot(3, 2, 1), mesh(g1l);
subplot(3, 2, 2), mesh(g2l);
conv_g1l_Icd = conv2(Icd, g1l, 'same');%这里用valid和same的结果不同
conv_g2l_Icd = conv2(Icd, g2l, 'same');
subplot(3, 2, 3), imagesc(conv_g1l_Icd);
subplot(3, 2, 4), imagesc(conv_g2l_Icd);
conv_g1l_Ia_grayd = conv2(Ia_grayd, g1l, 'same');
conv_g2l_Ia_grayd = conv2(Ia_grayd, g2l, 'same');
subplot(3, 2, 5), imagesc(conv_g1l_Ia_grayd);
subplot(3, 2, 6), imagesc(conv_g2l_Ia_grayd);
Q6 = conv_g1l_Icd(22, 41:45);
%% Q7: Multi-Scale Representations: The Gaussian Image Pyramid
gQ7 = fspecial('gaussian', 9, 1.5);

Temp = Ia_grayd;
figure
for i = 1:3
    conv_gQ7_Ia_grayd = conv2(Temp, gQ7, 'same');
    subplot(2, 2, i), imagesc(conv_gQ7_Ia_grayd)
    Temp = imresize(conv_gQ7_Ia_grayd, 0.5, 'nearest');
end
Q71 = Temp(20, 11);
Q72 = Temp(15, 20);
%% Q8: Multi-Scale Representations: The Laplacian Image Pyramid
g = fspecial('gaussian', 9, 1.5);

pyamid = Ia_grayd;
figure, clf
for i = 1:4
    conv_pyamid = conv2(pyamid, g, 'same');
    result = pyamid - conv_pyamid;
    subplot(2, 2, i), imagesc(result), colorbar;
    pyamid = imresize(conv_pyamid, 0.5, 'nearest');
end
Q81 = pyamid(5, 12);
Q82 = pyamid(11, 10);