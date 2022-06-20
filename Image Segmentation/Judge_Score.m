function [ score ] = Judge_Score( Ia, seg )
% Ia是人工描绘轮廓，seg是分割图像
Ia = im2bw(Ia, 0.5);
[m_Ia, n_Ia] = size(Ia);
[m_seg, n_seg] = size(seg);
% %% 随机选取118个点计算
% N = 118;% 共检测108个像素
% rand_row = ceil(rand(1, N)*m_Ia);
% rand_col = ceil(rand(1, N)*n_Ia);
% score = 0;
% for i = 1:N
%     score = score + sum(sum(abs(Ia(rand_row(i), rand_col(i))  ...
%     - seg(rand_row(i), rand_col(i))))) / N;
% end
% %% 全图像计算
% score = sum(sum(abs(Ia - seg))) ;% / (m_Ia * n_Ia);
%% f1 score = 2TP / (2TP + FP + FN)
% TP、True Positive   真阳性：预测为正，实际也为正
% FP、False Positive  假阳性：预测为正，实际为负
% FN、False Negative 假阴性：预测与负、实际为正
TP = 0; FP = 0; FN = 0; 
for i = 1:m_Ia
    for j = 1:n_Ia
        if seg(i, j) == 1 && Ia(i, j) == 1
            TP = TP + 1;
        end
        if seg(i, j) == 1 && Ia(i, j) == 0
            FP = FP + 1;
        end
        if seg(i, j) == 0 && Ia(i, j) == 1
            FN = FN + 1;
        end
    end
end
score = 2*TP / (2*TP + FP + FN);
end

