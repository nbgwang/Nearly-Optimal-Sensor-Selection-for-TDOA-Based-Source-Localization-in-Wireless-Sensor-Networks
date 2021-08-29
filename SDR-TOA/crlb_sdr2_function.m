%��ע�������и������������������޹�
%������������²���ʹ��
%����TOAת���ķ���������2��
%-------------------------��������б�------------------------
%  sensor_number������������
%  Noise_variance�����������������ķ�������м��������������ͬ��
%  sensor_sel_number��ѡ��Ĵ���������
%  Source_coordinate���ƶ�Ŀ��Դ�ĵ�ǰ���꣨��������
%  SENSOR_coordinate�������������꣨��������
%  dim:ά��

%-------------------------���ز����б�------------------------
%  crlb_01s��δ��˹�������CRLB
%  crlb_radoms����˹��������CRLB
%  cvx_cputime��cvx����ʱ��

%-------------------------�����б�------------------------
%  NOISE_covariance_matrix������������������Э�������
%  Select:������ѡ������
%  a:��С����ֵ��ϵ��������������Ӱ�죩
%----------------------------------------------------------
%��������
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
% %͹�Ż�
%������ѡ��������˹�����
T=sdr2cvx_unsensor(MEASUREMENT_matrix,NOISE_covariance_matrix,sensor_sel_number,sensor_number,dim);
crlb_GR=T(2);
crlb_unGR=T(1);
%---------------------------