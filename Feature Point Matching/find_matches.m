% ʹ�÷���
% ��ȡͼ��I1=im2double(imread('graffiti/img1.png')); I2=im2double(imread('graffiti/img2.png'));
% I1�е�Ŀ�������㼰λ�ã�detector = detectHarrisFeatures(rgb2gray(I1)); pos1=detector.Location;
% ������ƥ�䣺 [ pos2, tform ] = find_matches( I1, pos1, I2 ) 
%% �������ã���I2���ҵ���I1��pos1λ��ƥ������ص�
function [ pos2, tform ] = find_matches( I1, pos1, I2 ) 
%% Step 1: ��ȡͼ����Ԥ����
% bark, bikes, graffiti, leuven, trees, wall
% I1 = im2double(imread('./training_images/graffiti/img1.png')); I2 = im2double(imread('./training_images/graffiti/img2.png'));
% ceiling, day_night, semper, venice
% I1 = im2double(imread('./training_images/day_night/img1.jpg')); I2 = im2double(imread('./training_images/day_night/img2.jpg'));
%% Step 2: ��ȡHarris������(Ҳ������ȡSURF������)
% detect Harris
I2Points = detectHarrisFeatures(rgb2gray(I2)); 
pos2 = I2Points.Location; 
I1Points = detectHarrisFeatures(rgb2gray(I1)); 
pos1_all = I1Points.Location; 
% detect SURF
I2Points = detectSURFFeatures(rgb2gray(I2)); 
pos2 = [pos2; I2Points.Location]; 
I1Points = detectSURFFeatures(rgb2gray(I1)); 
pos1_all = [pos1_all; I1Points.Location]; 
% % detect BRISK
I2Points = detectBRISKFeatures(rgb2gray(I2)); 
pos2 = [pos2; I2Points.Location]; 
I1Points = detectBRISKFeatures(rgb2gray(I1)); 
pos1_all = [pos1_all; I1Points.Location]; 
% % detect FAST
% I2Points = detectFASTFeatures(rgb2gray(I2)); 
% pos2 = [pos2; I2Points.Location]; 
% I1Points = detectFASTFeatures(rgb2gray(I1)); 
% pos1_all = [pos1_all; I1Points.Location];
figure(3)
subplot(1, 2, 1)
imshow(I1);
hold on;
plot(pos1(:,1), pos1(:,2), 'go');
subplot(1, 2, 2)
imshow(I2);
hold on;
plot(pos2(:,1), pos2(:,2), 'ro');


% % detect MSER
% I2Points = detectMSERFeatures(rgb2gray(I2)); 
% pos2 = [pos2; I2Points.Location]; 
% % detect MinEigen
% I2Points = detectMinEigenFeatures(rgb2gray(I2)); 
% pos2 = [pos2; I2Points.Location]; 
%% Step 3: ��������������ͼ�����������(64ά��SURF����������64ά����128ά,��ƥ���ʱ��Ҫ�����п��ܼ������ң������ᱣ֤�ҵ��㹻��ĵ����ں���ƥ��)
[I1Features, I1Points] = extractFeatures(rgb2gray(I1), pos1_all, 'Method', 'SURF');
[I2Features, I2Points] = extractFeatures(rgb2gray(I2), pos2, 'Method', 'SURF');%, 'Method', 'SURF'
%% Step 4: ����������ƥ��ԣ���outliers��
Pairs = matchFeatures(I1Features, I2Features, 'MatchThreshold' , 50); %���ﻹ�и�����'MatchThreshold'����ֵƥ���֮��ľ���
matchedI1Points = I1Points(Pairs(:, 1), :);
matchedI2Points = I2Points(Pairs(:, 2), :);
% figure(1);
% showMatchedFeatures(I1, I2, matchedI1Points,matchedI2Points, 'montage');
% title('Putatively Matched Points (Including Outliers)'); 
% pause(1)
%% Step 5: Ԥ�����任��ȥ��������仯��outliers
[tform,inlierI1Points, inlierI2Points]  = estimateGeometricTransform(matchedI1Points, matchedI2Points, 'projective', 'MaxDistance', 5);
figure(1);
showMatchedFeatures(I1, I2, inlierI1Points, inlierI2Points, 'montage');
title('Matched Points (Inliers Only)');
pause(1)
%% Step 6: �õõ��ķ���任���󣬼���pos1�任���λ��pos2
pos2 = transformPointsForward(tform, pos1);
figure(2);
subplot(1, 2, 1)
imshow(I1);
hold on;
plot(pos1(: ,1), pos1(:,2),'+y');title('Img1 interest points(pos1)')
hold on
plot(pos1(119,1), pos1(119,2),'*g');
subplot(1, 2, 2)
imshow(I2);
hold on;
plot(pos2(: ,1), pos2(:,2),'+y');title('Img2 interest points(pos2)')
hold on
plot(pos2(119,1), pos2(119,2),'*r');
end

