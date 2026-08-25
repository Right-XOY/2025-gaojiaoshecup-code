function occ = q2_occlusion(M, C, R, params)
% q2_occlusion  判断导弹M处的烟幕云团是否将圆柱真目标完全遮蔽
%
% 判定准则（利用烟幕球体的凸性简化，见问题1结论）：
%   只需判断圆柱体上、下两个底面圆周是否被完全遮蔽，
%   即可判定整个圆柱体是否被完全遮蔽。
%
%   其中"圆周上一点P被遮蔽"定义为：
%     云团球（中心C、半径R）与线段 MP 相交，即
%     C 到线段 MP 的最近距离 <= R（垂足参数截断到 [0,1]，
%     允许云团包裹导弹端点或目标端点的情况）。
%
% 输入：
%   M      - 导弹位置，3x1 列向量 (m)
%   C      - 烟幕云团中心，3x1 列向量 (m)
%   R      - 烟幕有效半径 (m)，本题取 10
%   params - 场景结构体（字段：cylCenter, cylR, cylTop, cylBottom, nPhi）
% 输出：
%   occ    - logical，true 表示圆柱体被完全遮蔽

% ---- 上、下底面圆周采样点（各 nPhi 个，覆盖整圆 360°）----
P_top = q2_circle_points(params.cylCenter, params.cylR, params.cylTop, params.nPhi);
P_bot = q2_circle_points(params.cylCenter, params.cylR, params.cylBottom, params.nPhi);

% ---- 两个圆周的所有采样点均被遮蔽，才判定为完全遮蔽 ----
occ = all(q2_points_occluded(M, P_top, C, R)) && ...
      all(q2_points_occluded(M, P_bot, C, R));
end

%% ====================== 子函数 ======================

function P = q2_circle_points(center, r, z, n)
% 生成以 (center(1), center(2)) 为轴心、半径为 r、高度为 z 的圆周上的 n 个采样点
phi = linspace(0, 2*pi, n+1);
phi(end) = [];                          % 去掉重复端点，保留 n 个均匀分布点
P = [center(1) + r*cos(phi);
     center(2) + r*sin(phi);
     z*ones(1, n)];                     % 3 x n 矩阵
end

function occ = q2_points_occluded(M, P, C, R)
% 对一组目标点 P (3xn)，逐个判断其视线是否被烟幕球遮挡
v  = P - M;                             % 导弹 -> 目标点方向 (3xn)
w  = C - M;                             % 导弹 -> 云团中心 (3x1)
dvv = sum(v .* v, 1);                   % |v|^2 (1xn)
wv  = w' * v;                           % w·v (1xn)
t   = wv ./ dvv;                        % 垂足参数：foot = M + t*v

% 到线段 MP 的最近距离：垂足参数截断到 [0,1]（允许云团包裹导弹/目标端点）
t    = max(0, min(1, t));
foot = M + t .* v;                      % 线段上最近点坐标 (3xn)
d    = sqrt(sum((C - foot).^2, 1));     % 云团中心到线段 MP 的最近距离 (1xn)

occ = d <= R;                           % 烟幕球与视线段相交即遮蔽
end
