clc
clear all
% 思路一：写多种算法（阈值，扩张，聚类，拟合），求与observers平均最小差异的算法，作为最优分割方法，输出
% 思路二：深度学习方法，这种方法效果肯定最好，但是需要大量时间训练网络
%（或者直接从github down别人的代码，写别人的思路，并加引用），只提交测试训练好的网络
% 思路三：改进聚类或者拟合方法中的一种，需要查文献找trick
% Opencv工具包
% 前景背景分离算法`1
% 分别对不同通道进行分割
% 通过sum(sum(abs(* - *)))跟四副轮廓图比较，差值最小的为最优算法
% 对不同通道进行高斯平滑之后再聚类
%% 总结
% 去除小洞函数：bwareaopen(seg,100,8)，参数依次为二值图像，小洞面积，‘不知道’
% 侵蚀膨胀桥梁函数：bwmorph(二值图像，方法)
% 滤波器函数：imfilter(灰度图像，算子，方法)
% 提取轮廓函数：edge(灰度图像，算子"canny"，阈值)
% 聚集函数：conv2
% 计算灰度阈值函数：graythresh(f)
% 基于MATLAB的RGB转YCBCR色彩空间转换
%% 主函数
for j = 1:12
I=imread(['training_images/im', num2str(j),'.jpg']);
seg = segment_image(I);
Score = 0;
Ia_num = 5;
for i = 1:Ia_num
    Ia = imread(['training_images/', 'im', num2str(j), 'seg', num2str(i), '.png']);
    Score = Score + Judge_Score(Ia, seg);
end
SCORE(j) = Score/Ia_num;
end
mean(SCORE)
