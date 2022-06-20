clc
clear
%% 读取图像
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('woods.png');
%% resize图像大小
Ib(401:end,401:end)=255; %用有补丁的图像
Ibsmall=imresize(Ib,0.5);
Iblarge=imresize(Ib,2);
figure(1), clf
subplot(2,2,1), imagesc(Ibsmall)
subplot(2,2,2), imagesc(Iblarge)
Iblarge=imresize(Ib,2,'bilinear');
Ibsmall=imresize(Ib,0.5,'nearest');
Ibsmall_167_98 = Ibsmall(167,98)
Ib_167_98 = Ib(167,98)
Iblarge_167_98 = Iblarge(167,98)
Ib_334_196 = Ib(334,196)
Iblarge_668_392 = Iblarge(668,392)

clc
clear
%% 读取图像
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('woods.png');
%% Question 2
Ib(401:end,401:end)=255; %题目要求用 modified to have a white rectangular region
figure(2), clf
imagesc(Ib); colormap('gray')
Ibd=im2double(Ib); %将图像转换成double类型
Ibdiffv=Ibd(1:end-1,:)-Ibd(2:end,:); %图像边缘（垂直）
figure(3), clf, imagesc(Ibdiffv);colormap('gray');colorbar
Ibdiffh=Ibd(:,1:end-1)-Ibd(:,2:end); %图像边缘（水平）
figure(4), clf, imagesc(Ibdiffh);colormap('gray');colorbar
Ibdiff=sqrt(Ibdiffh(1:end-1,:).^2+Ibdiffv(:,1:end-1).^2); %行列边缘整合
figure(5), clf, imagesc(Ibdiff); colormap('gray'); colorbar
bw=im2bw(Ibdiff,0.075); %转换成二元图像
figure(6), clf, imagesc(bw); colormap('gray'); colorbar
Ibdiffv_473_400 = Ibdiffv(473, 400)
Ibdiffh_473_400 = Ibdiffh(473, 400)
Ibdiff_473_400 = Ibdiff(473, 400)
bw_473_400 = bw(473, 400)

clc
clear
%% 读取图像
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('woods.png');
%% question 3 图像的冗余，相关系数
Ib(401:end,401:end)=255;
offsets = [0:30];
for offset = offsets
    simb(offset+1) = corr2(Ib(1:end - offset,:), Ib(1+offset:end, :));
    simc(offset+1) = corr2(Ic(1:end - offset,:), Ic(1+offset:end, :));
end
figure(7),clf
plot(offsets, simb);
simb_5 = simb(5);
simb_27 = simb(27);
hold on
plot(offsets, simc);
simc_5 = simc(5);
simc_27 = simc(27);
legend({'elephant', 'woods'});
xlabel('shift');
ylabel('correlation coefficient')

clc
clear
%% 读取图像
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('woods.png');
%% question 4 为了模拟视网膜神经节细胞对图像所有不同部分的响应，将高斯差分模板与图像进行卷积
Ib(401:end,401:end)=255; %不知道这里要不要采用补丁的大象图像
Ibd=im2double(Ib); %将图像转换成double类型
Icd=im2double(Ic); %将图像转换成double类型
dog=fspecial('gaussian',9,1)-fspecial('gaussian',9,1.5);
Ibdog=conv2(Ibd,dog,'same');
Icdog=conv2(Icd,dog,'same');
% figure(8), clf, imagesc(Ibdog); colormap('gray'); colorbar
% figure(8), clf, imagesc(Icdog); colormap('gray'); colorbar
offsets = [0:30];
for offset = offsets
    simb2(offset+1) = corr2(Ibdog(1:end - offset,:), Ibdog(1+offset:end, :));
    simc2(offset+1) = corr2(Icdog(1:end - offset,:), Icdog(1+offset:end, :));
end
figure(7),clf
plot(offsets, simb2);
simb2_9 = simb2(9); % [0.318963406598446]
simb2_22 = simb2(22); % [0.311735710100967]
hold on
plot(offsets, simc2);
simc2_9 = simc2(9); % [0.0635676992899844]
simc2_22 = simc2(22); % [0.104447645755953]
legend({'elephant', 'woods'});
xlabel('shift');
ylabel('correlation coefficient')

clc
clear
%% 读取图像
Ia=imread('rooster.jpg');
Ib=imread('elephant.png');
Ic=imread('woods.png');
%% Question 5 颜色检测
Iad=im2double(Ia);
gc=fspecial('gaussian',9,1); %中心
gs=fspecial('gaussian',9,1.5); %surround
IaRG=conv2(Iad(:,:,1),gc,'same')-conv2(Iad(:,:,2),gs,'same'); %红开绿关
IaGR=conv2(Iad(:,:,2),gc,'same')-conv2(Iad(:,:,1),gs,'same'); %绿开红关
IaBY=conv2(Iad(:,:,3),gc,'same')-conv2(mean(Iad(:,:,1:2),3),gs,'same'); %蓝开黄关
IaYB=conv2(mean(Iad(:,:,1:2),3),gc,'same')-conv2(Iad(:,:,3),gs,'same'); %黄开蓝关
IaRG_244_385 = IaRG(244, 385)
IaGR_244_385 = IaGR(244, 385)
IaBY_244_385 = IaBY(244, 385)
IaYB_244_385 = IaYB(244, 385)