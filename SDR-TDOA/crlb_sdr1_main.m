clc;clear;close all;
sensor_number=20;
sensor_sel_number=5;
dim=2;
Source_coordinate=unifrnd(0,sensor_number*2.5,dim,1);
SENSOR_coordinate=unifrnd(0,sensor_number*2.5,dim,sensor_number);
for k=1:sensor_number
    Range(k,:)=norm((Source_coordinate-SENSOR_coordinate(:,k)),2);
end
RANGE=Range*ones(1,dim);
MEASUREMENT_matrix=((Source_coordinate*ones(1,sensor_number))'-(SENSOR_coordinate)')./RANGE;
%-------------------------------------
%����
Q=eye(sensor_number);
a=0.1;
Q0=Q-a*eye(sensor_number);
%-------------------------------------
%��������
C=MEASUREMENT_matrix'/(Q0)*MEASUREMENT_matrix;
B=Q0\MEASUREMENT_matrix;
%-------------------------------------
%͹�Ż�
%ѡ�񴫸�����CRLB
CRLB_tdoa=sdr1cvx_unsensor(Q0,B,C,a,sensor_sel_number,sensor_number,dim);
%δ������˹���������
crlb_sdr1uGR=CRLB_tdoa(1);
%������˹���������
crlb_sdr1GR=CRLB_tdoa(2);