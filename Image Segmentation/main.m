clc
clear all
% ˼·һ��д�����㷨����ֵ�����ţ����࣬��ϣ�������observersƽ����С������㷨����Ϊ���ŷָ�������
% ˼·�������ѧϰ���������ַ���Ч���϶���ã�������Ҫ����ʱ��ѵ������
%������ֱ�Ӵ�github down���˵Ĵ��룬д���˵�˼·���������ã���ֻ�ύ����ѵ���õ�����
% ˼·�����Ľ����������Ϸ����е�һ�֣���Ҫ��������trick
% Opencv���߰�
% ǰ�����������㷨`1
% �ֱ�Բ�ͬͨ�����зָ�
% ͨ��sum(sum(abs(* - *)))���ĸ�����ͼ�Ƚϣ���ֵ��С��Ϊ�����㷨
% �Բ�ͬͨ�����и�˹ƽ��֮���پ���
%% �ܽ�
% ȥ��С��������bwareaopen(seg,100,8)����������Ϊ��ֵͼ��С�����������֪����
% ��ʴ��������������bwmorph(��ֵͼ�񣬷���)
% �˲���������imfilter(�Ҷ�ͼ�����ӣ�����)
% ��ȡ����������edge(�Ҷ�ͼ������"canny"����ֵ)
% �ۼ�������conv2
% ����Ҷ���ֵ������graythresh(f)
% ����MATLAB��RGBתYCBCRɫ�ʿռ�ת��
%% ������
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
