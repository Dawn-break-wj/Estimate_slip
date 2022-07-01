%%
clc,clear;
%% Load Data
load('sliding5.mat')

%% 在采集的数据中提取所需要的信息
ay=sliding5.Y(1).Data;%采集的纵向加速度
Current_Speed=sliding5.Y(17).Data;%车速
LR_Wheel_Speed=sliding5.Y(24).Data/36;%左轮速
RR_Wheel_Speed=sliding5.Y(27).Data/36;%右轮速
Motor_Torque=(sliding5.Y(26).Data+18000);
T=0.02;

%% 左轮速度数据处理（Kalman Filter）
%kalman滤波（二维，x,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化参数
N=length(ay);
Q=[0.05 0;0 0.05];
R=1;
w=sqrt(Q)*randn(2,N);
v=sqrt(R)*randn(1,N);
A=[1 T;0 1];%状态转移矩阵
B=[0;0];%控制矩阵
U=0;%控制变量
H=[0,1];%观测矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%分配空间，x,p0,z,xkf
x=zeros(2,N);%物体真实状态
x(:,1)=[0;0];%初始位移和速度
p0=[0 0;0 0];%初始误差
z=zeros(1,N);%观测值
z(1)=LR_Wheel_Speed(:,1);%观测真实值,第一列的列向量，
xkf=zeros(2,N);%kalman估计状态
xkf(:,1)=x(:,1);%kalman估计状态初始化,第一列的列向量，
err_p=zeros(N,2);
err_p(1,1)=p0(1,1);
err_p(1,2)=p0(2,2);
I=eye(2);%2*2单位矩阵表明二维系统
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kalman迭代过程
for k=2:N
    x(:,k)=A*x(:,k-1)+B*U+w(k);%模型
    z(k)=LR_Wheel_Speed(:,k);%观测值
    x_pre=A*xkf(:,k-1)+B*U;%①
    p_pre=A*p0*A'+Q;%②
    kg=p_pre*H'/(H*p_pre*H'+R);%③
    xkf(:,k)=x_pre+kg*(z(k)-H*x_pre);%④
    p0=(I-kg*H)*p_pre;%⑤
    err_p(k,1)=p0(1,1);%位移误差均方值
    err_p(k,2)=p0(2,2);%速度误差均方值
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%误差计算
messure_err_v=zeros(1,N);
kalman_err_x=zeros(1,N);
kalman_err_v=zeros(1,N);
for k=1:N
    messure_err_v(k)=z(k)-v(1,k);%速度的测量误差
    kalman_err_x(k)=xkf(1,k)-x(1,k);%kalman估计位移与真实位移之间偏差
    kalman_err_v(k)=xkf(2,k)-x(2,k);%kalman估计速度与真实速度之间偏差
end
% 绘图
figure;
plot(xkf(2,:));
hold on
plot(Current_Speed);
hold on
plot(LR_Wheel_Speed)
LR_Wheel_Speed=xkf(2,:);
hold on
%% 左轮车速的打滑率计算
[size_ax,size_ay]=size(ay);
Size_Data=max(size_ax,size_ay);
% 打滑状态判断和参考轮速计算
for i=1:Size_Data-1
    if (LR_Wheel_Speed(1,i)*1.5<LR_Wheel_Speed(1,i+1) && ay(1,i+1)>0) ||...
            (LR_Wheel_Speed(1,i)>1.5*LR_Wheel_Speed(1,i+1) && ay(1,i+1)<0)
        LR_Wheel_Speed(1,i+1)=LR_Wheel_Speed(1,i)+ay(1,i)*T;
    end
end
    
for k=1:Size_Data
    s_LR(1,k)=(abs(Current_Speed(1,k)-LR_Wheel_Speed(1,k)))/max(Current_Speed(1,k),LR_Wheel_Speed(1,k));
end

%% 右轮速度数据处理（Kalman Filter）
%kalman滤波（二维，x,v,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化参数
N=length(ay);
Q=[0.05 0 0;0 0.05 0;0 0 0.05];
R=1;
w=sqrt(Q)*randn(3,N);
v=sqrt(R)*randn(1,N);
A=[1 T T*T/2;0 1 T;0 0 1];%状态转移矩阵
B=[0;0;0];%控制矩阵
U=0;%控制变量
H=[0,1,1];%观测矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%分配空间，x,p0,z,xkf
x=zeros(3,N);%物体真实状态
x(:,1)=[0;0;0];%初始位移和速度
p0=[0 0 0;0 0 0;0 0 0];%初始误差
z=zeros(1,N);%观测值
z(1)=RR_Wheel_Speed(:,1);%观测真实值,第一列的列向量，
xkf=zeros(3,N);%kalman估计状态
xkf(:,1)=x(:,1);%kalman估计状态初始化,第一列的列向量，
err_p=zeros(N,2);
err_p(1,1)=p0(1,1);
err_p(1,2)=p0(2,2);
I=eye(3);%2*2单位矩阵表明二维系统
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kalman迭代过程
for k=2:N
    x(:,k)=A*x(:,k-1)+B*U+w(k);%模型
    z(k)=LR_Wheel_Speed(:,k);%观测值
    x_pre=A*xkf(:,k-1)+B*U;%①
    p_pre=A*p0*A'+Q;%②
    kg=p_pre*H'/(H*p_pre*H'+R);%③
    xkf(:,k)=x_pre+kg*(z(k)-H*x_pre);%④
    p0=(I-kg*H)*p_pre;%⑤
    err_p(k,1)=p0(1,1);%位移误差均方值
    err_p(k,2)=p0(2,2);%速度误差均方值
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%误差计算
messure_err_x=zeros(1,N);
kalman_err_x=zeros(1,N);
kalman_err_v=zeros(1,N);
for k=1:N
    messure_err_v(k)=z(k)-v(1,k);%速度的测量误差
    kalman_err_x(k)=xkf(1,k)-x(1,k);%kalman估计位移与真实位移之间偏差
    kalman_err_v(k)=xkf(2,k)-x(2,k);%kalman估计速度与真实速度之间偏差
end

RR_Wheel_Speed=xkf(2,:);

%% 右轮车速的打滑率计算
% 打滑状态判断和参考轮速计算
for i=1:Size_Data-1
    if (RR_Wheel_Speed(1,i)*1.5<RR_Wheel_Speed(1,i+1) && ay(1,i+1)>0) ||...
            (RR_Wheel_Speed(1,i)>1.5*RR_Wheel_Speed(1,i+1) && ay(1,i+1)<0)
        RR_Wheel_Speed(1,i+1)=RR_Wheel_Speed(1,i)+ay(1,i)*T;
    end
end
    
for k=1:Size_Data
    s_RR(1,k)=(abs(Current_Speed(1,k)-RR_Wheel_Speed(1,k)))/max(Current_Speed(1,k),RR_Wheel_Speed(1,k));
end
%% 计算结果处理
%  如果车速小于0.1m/s或者轮速小于0.1m/s，直接令s=0
for j = 1:Size_Data
    if Current_Speed(j)<0.1 || LR_Wheel_Speed(j)<0.1 || RR_Wheel_Speed(j)<0.1
        s_LR(j)= 0;
        s_RR(j)= 0;
    end
end

%% 利用递归最小二乘法(RLS)计算滑移斜率
Wheel_Torque = Motor_Torque*12/2;
F=Wheel_Torque/0.247;%车轮滚动半径为247mm

% 根据F=ma,地面作用在车轮上的牵引力-ma=摩擦斜率*滑移率*mg=mu*mg
mu = abs((2*F - 240*ay))/(240*9.8);
theta=(4*F)/(240*9.8);

%% RLS计算
%RLS引入数据
phi = s_LR(1,:)';%设置测量回归矢量phi
b = mu(1,:)';%测量输出量
x = 0;%未知参数（即要求的数）
I = eye(1,1);
P = (10^6) * I;%协方差矩阵
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