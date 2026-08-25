function occ = q3_occlusion_window(Mw, Cw, R, params)
% q3_occlusion_window  对一段连续时间窗，向量化判定圆柱真目标是否被烟幕云团完全遮蔽
%
% 判定准则（与问题2一致，利用凸性简化）：
%   圆柱上、下两个底面圆周的所有采样点均被烟幕球遮挡，
%   即判定整个圆柱体被完全遮蔽。
%   其中"点 P 被遮挡" = 烟幕球（中心 C、半径 R）与线段 MP 相交，
%   即 C 到线段 MP 的最近距离 <= R（垂足参数截断到 [0,1]，
%   允许云团包裹导弹端点或目标端点的情况）。
%
% 输入：
%   Mw    - 导弹位置序列，3 x T (m)
%   Cw    - 烟幕云团中心序列，3 x T (m)，与 Mw 时间一一对应
%   R     - 烟幕有效半径 (m)
%   params- 场景结构体（cylCenter, cylR, cylTop, cylBottom, nPhi）
% 输出：
%   occ   - 1 x T 逻辑向量，occ(k)=true 表示第 k 个时刻完全遮蔽

nPhi = params.nPhi;
phi = linspace(0, 2*pi, nPhi + 1);
phi(end) = [];                                   % nPhi 个均匀采样角
P_top = [params.cylCenter(1) + params.cylR*cos(phi);
         params.cylCenter(2) + params.cylR*sin(phi);
         params.cylTop*ones(1, nPhi)];           % 上底面圆周 (3 x nPhi)
P_bot = [params.cylCenter(1) + params.cylR*cos(phi);
         params.cylCenter(2) + params.cylR*sin(phi);
         params.cylBottom*ones(1, nPhi)];        % 下底面圆周

T = size(Mw, 2);
% 扩展为 3 x 1 x T，便于隐式扩展与圆周采样点做叉乘
M3 = reshape(Mw, 3, 1, T);
C3 = reshape(Cw, 3, 1, T);

occ = true(1, 1, T);                     % 1 x 1 x T，与 all(...,2) 结果同形状，
                                         % 避免 1 x T 与 1 x 1 x T 按 & 隐式扩展成 1 x T x T
for P = {P_top, P_bot}
    Pn = reshape(P{1}, 3, [], 1);                % 3 x nPhi x 1
    d = Pn - M3;                                 % 导弹->采样点方向 (3 x nPhi x T)
    w = C3 - M3;                                 % 导弹->云团中心 (3 x 1 x T)
    dvv = sum(d.^2, 1);                          % (1 x nPhi x T)
    wv  = sum(w .* d, 1);                        % (1 x nPhi x T)
    tp  = wv ./ dvv;                             % 垂足参数
    tp  = max(0, min(1, tp));                    % 截断到 [0,1]（到线段距离，含端点）
    foot = M3 + tp .* d;                         % 线段上最近点 (3 x nPhi x T)
    dist = sqrt(sum((C3 - foot).^2, 1));         % C 到线段 MP 最近距离 (1 x nPhi x T)
    occ = occ & all(dist <= R, 2);               % 该圆上所有点均被遮挡
end

occ = squeeze(occ);                              % 1 x T（T=1 时为标量）
if size(occ, 1) > 1, occ = occ'; end
if ~isrow(occ) && numel(occ) ~= 1                % 形状保护：必须为行向量或标量
    error('q3_occlusion_window 输出形状异常 %s，请 clear functions 后重试', ...
          mat2str(size(occ)));
end
end
