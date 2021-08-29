%备注：仿真中各传感器测量噪声均无关
%测量噪声相关下不能使用
%基于TOA转换的方法（方法2）
%-------------------------输入参数列表------------------------
%  sensor_number：传感器总数
%  Noise_variance：传感器测量噪声的方差（仿真中假设噪声方差均相同）
%  sensor_sel_number：选择的传感器数量
%  Source_coordinate：移动目标源的当前坐标（列向量）
%  SENSOR_coordinate：传感器的坐标（列向量）
%  dim:维数

%-------------------------返回参数列表------------------------
%  crlb_01s：未高斯随机化的CRLB
%  crlb_radoms：高斯随机化后的CRLB
%  cvx_cputime：cvx运行时间

%-------------------------参数列表------------------------
%  NOISE_covariance_matrix：传感器测量噪声的协方差矩阵
%  Select:传感器选择向量
%  a:最小特征值的系数（本程序中无影响）
%----------------------------------------------------------
%数据生成
warning off
clc;clear;close all;
sensor_number=20;
sensor_sel_number=5;
dim=2;
Source_coordinate=unifrnd(0,sensor_number*2.5,dim,1);
SENSOR_coordinate=unifrnd(0,sensor_number*2.5,dim,sensor_number);
NOISE_covariance_matrix=eye(sensor_number);
for k=1:sensor_number
    Range(k,:)=norm((Source_coordinate-SENSOR_coordinate(:,k)),2);
end
RANGE=Range*ones(1,dim);
MEASUREMENT_matrix=[((Source_coordinate*ones(1,sensor_number))'-(SENSOR_coordinate)')./RANGE,ones(sensor_number,1)];
%---------------------------
% %凸优化
%传感器选择向量高斯随机化
T=sdr2cvx_unsensor(MEASUREMENT_matrix,NOISE_covariance_matrix,sensor_sel_number,sensor_number,dim);
crlb_GR=T(2);
crlb_unGR=T(1);
%---------------------------