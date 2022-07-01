%%
clc,clear;
%% Load Data
load('sliding5.mat')

%% �ڲɼ�����������ȡ����Ҫ����Ϣ
ay=sliding5.Y(1).Data;%�ɼ���������ٶ�
Current_Speed=sliding5.Y(17).Data;%����
LR_Wheel_Speed=sliding5.Y(24).Data/36;%������
RR_Wheel_Speed=sliding5.Y(27).Data/36;%������
Motor_Torque=(sliding5.Y(26).Data+18000);
T=0.02;

%% �����ٶ����ݴ���Kalman Filter��
%kalman�˲�����ά��x,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ������
N=length(ay);
Q=[0.05 0;0 0.05];
R=1;
w=sqrt(Q)*randn(2,N);
v=sqrt(R)*randn(1,N);
A=[1 T;0 1];%״̬ת�ƾ���
B=[0;0];%���ƾ���
U=0;%���Ʊ���
H=[0,1];%�۲����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ռ䣬x,p0,z,xkf
x=zeros(2,N);%������ʵ״̬
x(:,1)=[0;0];%��ʼλ�ƺ��ٶ�
p0=[0 0;0 0];%��ʼ���
z=zeros(1,N);%�۲�ֵ
z(1)=LR_Wheel_Speed(:,1);%�۲���ʵֵ,��һ�е���������
xkf=zeros(2,N);%kalman����״̬
xkf(:,1)=x(:,1);%kalman����״̬��ʼ��,��һ�е���������
err_p=zeros(N,2);
err_p(1,1)=p0(1,1);
err_p(1,2)=p0(2,2);
I=eye(2);%2*2��λ���������άϵͳ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kalman��������
for k=2:N
    x(:,k)=A*x(:,k-1)+B*U+w(k);%ģ��
    z(k)=LR_Wheel_Speed(:,k);%�۲�ֵ
    x_pre=A*xkf(:,k-1)+B*U;%��
    p_pre=A*p0*A'+Q;%��
    kg=p_pre*H'/(H*p_pre*H'+R);%��
    xkf(:,k)=x_pre+kg*(z(k)-H*x_pre);%��
    p0=(I-kg*H)*p_pre;%��
    err_p(k,1)=p0(1,1);%λ��������ֵ
    err_p(k,2)=p0(2,2);%�ٶ�������ֵ
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������
messure_err_v=zeros(1,N);
kalman_err_x=zeros(1,N);
kalman_err_v=zeros(1,N);
for k=1:N
    messure_err_v(k)=z(k)-v(1,k);%�ٶȵĲ������
    kalman_err_x(k)=xkf(1,k)-x(1,k);%kalman����λ������ʵλ��֮��ƫ��
    kalman_err_v(k)=xkf(2,k)-x(2,k);%kalman�����ٶ�����ʵ�ٶ�֮��ƫ��
end
% ��ͼ
figure;
plot(xkf(2,:));
hold on
plot(Current_Speed);
hold on
plot(LR_Wheel_Speed)
LR_Wheel_Speed=xkf(2,:);
hold on
%% ���ֳ��ٵĴ��ʼ���
[size_ax,size_ay]=size(ay);
Size_Data=max(size_ax,size_ay);
% ��״̬�жϺͲο����ټ���
for i=1:Size_Data-1
    if (LR_Wheel_Speed(1,i)*1.5<LR_Wheel_Speed(1,i+1) && ay(1,i+1)>0) ||...
            (LR_Wheel_Speed(1,i)>1.5*LR_Wheel_Speed(1,i+1) && ay(1,i+1)<0)
        LR_Wheel_Speed(1,i+1)=LR_Wheel_Speed(1,i)+ay(1,i)*T;
    end
end
    
for k=1:Size_Data
    s_LR(1,k)=(abs(Current_Speed(1,k)-LR_Wheel_Speed(1,k)))/max(Current_Speed(1,k),LR_Wheel_Speed(1,k));
end

%% �����ٶ����ݴ���Kalman Filter��
%kalman�˲�����ά��x,v,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ������
N=length(ay);
Q=[0.05 0 0;0 0.05 0;0 0 0.05];
R=1;
w=sqrt(Q)*randn(3,N);
v=sqrt(R)*randn(1,N);
A=[1 T T*T/2;0 1 T;0 0 1];%״̬ת�ƾ���
B=[0;0;0];%���ƾ���
U=0;%���Ʊ���
H=[0,1,1];%�۲����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ռ䣬x,p0,z,xkf
x=zeros(3,N);%������ʵ״̬
x(:,1)=[0;0;0];%��ʼλ�ƺ��ٶ�
p0=[0 0 0;0 0 0;0 0 0];%��ʼ���
z=zeros(1,N);%�۲�ֵ
z(1)=RR_Wheel_Speed(:,1);%�۲���ʵֵ,��һ�е���������
xkf=zeros(3,N);%kalman����״̬
xkf(:,1)=x(:,1);%kalman����״̬��ʼ��,��һ�е���������
err_p=zeros(N,2);
err_p(1,1)=p0(1,1);
err_p(1,2)=p0(2,2);
I=eye(3);%2*2��λ���������άϵͳ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kalman��������
for k=2:N
    x(:,k)=A*x(:,k-1)+B*U+w(k);%ģ��
    z(k)=LR_Wheel_Speed(:,k);%�۲�ֵ
    x_pre=A*xkf(:,k-1)+B*U;%��
    p_pre=A*p0*A'+Q;%��
    kg=p_pre*H'/(H*p_pre*H'+R);%��
    xkf(:,k)=x_pre+kg*(z(k)-H*x_pre);%��
    p0=(I-kg*H)*p_pre;%��
    err_p(k,1)=p0(1,1);%λ��������ֵ
    err_p(k,2)=p0(2,2);%�ٶ�������ֵ
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������
messure_err_x=zeros(1,N);
kalman_err_x=zeros(1,N);
kalman_err_v=zeros(1,N);
for k=1:N
    messure_err_v(k)=z(k)-v(1,k);%�ٶȵĲ������
    kalman_err_x(k)=xkf(1,k)-x(1,k);%kalman����λ������ʵλ��֮��ƫ��
    kalman_err_v(k)=xkf(2,k)-x(2,k);%kalman�����ٶ�����ʵ�ٶ�֮��ƫ��
end

RR_Wheel_Speed=xkf(2,:);

%% ���ֳ��ٵĴ��ʼ���
% ��״̬�жϺͲο����ټ���
for i=1:Size_Data-1
    if (RR_Wheel_Speed(1,i)*1.5<RR_Wheel_Speed(1,i+1) && ay(1,i+1)>0) ||...
            (RR_Wheel_Speed(1,i)>1.5*RR_Wheel_Speed(1,i+1) && ay(1,i+1)<0)
        RR_Wheel_Speed(1,i+1)=RR_Wheel_Speed(1,i)+ay(1,i)*T;
    end
end
    
for k=1:Size_Data
    s_RR(1,k)=(abs(Current_Speed(1,k)-RR_Wheel_Speed(1,k)))/max(Current_Speed(1,k),RR_Wheel_Speed(1,k));
end
%% ����������
%  �������С��0.1m/s��������С��0.1m/s��ֱ����s=0
for j = 1:Size_Data
    if Current_Speed(j)<0.1 || LR_Wheel_Speed(j)<0.1 || RR_Wheel_Speed(j)<0.1
        s_LR(j)= 0;
        s_RR(j)= 0;
    end
end

%% ���õݹ���С���˷�(RLS)���㻬��б��
Wheel_Torque = Motor_Torque*12/2;
F=Wheel_Torque/0.247;%���ֹ����뾶Ϊ247mm

% ����F=ma,���������ڳ����ϵ�ǣ����-ma=Ħ��б��*������*mg=mu*mg
mu = abs((2*F - 240*ay))/(240*9.8);
theta=(4*F)/(240*9.8);

%% RLS����
%RLS��������
phi = s_LR(1,:)';%���ò����ع�ʸ��phi
b = mu(1,:)';%���������
x = 0;%δ֪��������Ҫ�������
I = eye(1,1);
P = (10^6) * I;%Э�������
kesi=0.98;
for k = 1:466
    phi_k = phi(k,:);
    Q1 = P*(phi_k');
    Q2 = kesi + phi_k * P * (phi_k');
    Q = Q1/Q2;
    x = x + Q * (b(k) - phi_k*x);
    P = 1/kesi*((I - Q*phi_k)*P);
    result2(:,k) = x;
    result1(k) = k;
end
result1 = result1';
%result = [result1; result2];
plot(result1, result2);
s_LR_10=s_LR*10;