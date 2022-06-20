function [seg] = segment_image(I)
%% �㷨1����Ե��ⷽ��
I = im2double(I);
    sigma = 0.01;
    freq = 0.01;
    aspect = 0.75;
    phase1 = 90;
    phase2 = 0;
    Idg = rgb2gray(I);
    Img = 0;
    for orientation = 0:10:180
        g1 = gabor2(sigma,freq,orientation,aspect,phase1);
        g2 = gabor2(sigma,freq,orientation,aspect,phase2);
        img1 = conv2(Idg, g1, 'same');
        img2 = conv2(Idg, g2, 'same');
        img = sqrt((img1.^2) + (img2.^2));
        Img = max(Img, img);
    end
    mask_mean = fspecial('average', [3 3]);
    Img = conv2(Img, mask_mean, 'same');
    seg = edge(Img, 'canny');%��ȡ��Ե
    seg = bwareaopen(seg, 100, 8) ;%ȥ��С��
    figure, imshow(seg), title('ȥ��С��');
end