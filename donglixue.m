    clc;clear;
%所有物理量都是国际单位
%mdh参数，不是dh参数
dh(1,:) = [    0, -0.285,      0,    pi, 0, 0];
dh(2,:) = [    0,      0,      0,  pi/2, 0, pi/2];
dh(3,:) = [    0,      0, -0.398,    pi, 0, 0];
dh(4,:) = [    0, -0.194,  -0.36,    pi, 0, -pi/2];
dh(5,:) = [    0,  -0.17,      0, -pi/2, 0, 0];
dh(6,:) = [    0,  -0.16,      0,  pi/2, 0, 0];

% 记第i个连杆由MDH建模方法决定的随体坐标系为{i}
% 记第i个连杆的质心坐标系为{ic}
% {ic}与{i}方向相同
Inertia = [[0.017, 0.466, 0.017, 0, -0.001, 0],[0, 0.043, 0.005], 5.675;
           [0.014, 0.124, 0.122, 0, 0, -0.009],[-0.298, 0, -0.126], 5.466;
           [0.01, 0.081, 0.076, 0, 0, -0.01],[-0.274, 0, 0.061],4.294;
           [0.005, 0.003, 0.005, 0, 0, 0],[0, -0.055, 0.001],2.732;
           [0.005, 0.003, 0.005, 0, 0, 0],[0, 0.055, 0.001],2.732;
           [0.0001, 0.0001, 0.0001, 0, 0, 0],[0, 0, 0.017], 0.28];      
%Inertia(i,:) = [[Ixx,Iyy,Izz,Ixy,Iyz,Ixz],[rx,ry,rz],m]
%[Ixx,Iyy,Izz,Ixy,Iyz,Ixz]:第i个连杆在{i}下的惯量
%[rx,ry,rz]：第i个连杆质心在{i}下的坐标
%m：连杆质量       
for i = 1:6
L(i) = Link(dh(i,:),'modified');
I = zeros(3,3); I(1,2:3) = Inertia(i,[4,6]);I(2,3) = Inertia(i,5);
I = I+I'; I = I + diag(Inertia(i,1:3)) ;
L(i).I = I;
L(i).r = Inertia(i,7:9); %L(i).r=[rx,ry,rz]
L(i).m = Inertia(i,10);
L(i).I =L(i).I + L(i).m*skew(L(i).r)^2; %平行轴定理,L(i).I是第i个连杆在{ic}下的惯量
L(i).Jm = 0;
end
robot = SerialLink(L,'name', 'qinzhe');

%取重力加速度g=9.81m/s2
g=9.81;

syms t1 t2 t3 t4 t5 t6 dt1 dt2 dt3 dt4 dt5 dt6 ddt1 ddt2 ddt3 ddt4 ddt5 ddt6;
T01=[1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t1) -sin(t1) 0 0;sin(t1) cos(t1) 0 0;0 0 1 -0.285;0 0 0 1];
T12=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[-sin(t2) -cos(t2) 0 0;cos(t2) -sin(t2) 0 0;0 0 1 0;0 0 0 1];
T23=[1 0 0 -0.398;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t3) -sin(t3) 0 0;sin(t3) cos(t3) 0 0;0 0 1 0;0 0 0 1];
T34=[1 0 0 -0.36;0 -1 0 0;0 0 -1 0;0 0 0 1]*[sin(t4) cos(t4) 0 0;-cos(t4) sin(t4) 0 0;0 0 1 -0.194;0 0 0 1];
T45=[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1]*[cos(t5) -sin(t5) 0 0;sin(t5) cos(t5) 0 0;0 0 1 -0.17;0 0 0 1];
T56=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[cos(t6) -sin(t6) 0 0;sin(t6) cos(t6) 0 0;0 0 1 -0.16;0 0 0 1];





T02=simplify(T01*T12);
T03=simplify(T02*T23);
T04=simplify(T03*T34);
T05=simplify(T04*T45);
T06=simplify(T05*T56);



%对于雅可比中的第i列，不应计算z(i-1)，而应是zi, 原因为mdh系与dh系不同，mdh系qi表征zi的转动,z0没有用
%mdh系中0系除了表示世界坐标系外没有其他用，没有7系
z=sym(zeros(3,6));
z(:,1)=T01(1:3,3);
z(:,2)=T02(1:3,3);
z(:,3)=T03(1:3,3);
z(:,4)=T04(1:3,3);
z(:,5)=T05(1:3,3);
z(:,6)=T06(1:3,3);  %修改后的z(:,i)表示zi

%修改的原理同上
P=sym(zeros(3,6));
P(:,1)=T01(1:3,4);
P(:,2)=T02(1:3,4);
P(:,3)=T03(1:3,4);
P(:,4)=T04(1:3,4);
P(:,5)=T05(1:3,4);
P(:,6)=T06(1:3,4);

%计算平动动能之时，应计算质心的速度而非坐标系原点的速度, 故雅可比中为cross(zi,Pc-P)而非cross(zi,P-P)
Pc1=T01*[L(1).r.';1];
Pc2=T02*[L(2).r.';1];
Pc3=T03*[L(3).r.';1];
Pc4=T04*[L(4).r.';1];
Pc5=T05*[L(5).r.';1];
Pc6=T06*[L(6).r.';1];


Jl1=sym(zeros(6,6));
for i=1:1
    Jl1(:,i)=[cross(z(:,i),Pc1(1:3,1)-P(:,i));z(:,i)];
end

Jl2=sym(zeros(6,6));
for i=1:2
    Jl2(:,i)=[cross(z(:,i),Pc2(1:3,1)-P(:,i));z(:,i)];
end


Jl3=sym(zeros(6,6));
for i=1:3
    Jl3(:,i)=[cross(z(:,i),Pc3(1:3,1)-P(:,i));z(:,i)];
end

Jl4=sym(zeros(6,6));
for i=1:4
    Jl4(:,i)=[cross(z(:,i),Pc4(1:3,1)-P(:,i));z(:,i)];
end


Jl5=sym(zeros(6,6));
for i=1:5
    Jl5(:,i)=[cross(z(:,i),Pc5(1:3,1)-P(:,i));z(:,i)];
end

Jl6=sym(zeros(6,6));
for i=1:6
    Jl6(:,i)=[cross(z(:,i),Pc6(1:3,1)-P(:,i));z(:,i)];
end



%动能和势能计算
q=[t1;t2;t3;t4;t5;t6];
q_dot=[dt1;dt2;dt3;dt4;dt5;dt6];
q_dotdot=[ddt1;ddt2;ddt3;ddt4;ddt5;ddt6];
Tk=sym(zeros(6,1));
Tk(1)=0.5*L(1).m*q_dot.'*Jl1(1:3,:).'*Jl1(1:3,:)*q_dot+0.5*q_dot.'*Jl1(4:6,:).'*T01(1:3,1:3)*L(1).I*T01(1:3,1:3).'*Jl1(4:6,:)*q_dot;
Tk(2)=0.5*L(2).m*q_dot.'*Jl2(1:3,:).'*Jl2(1:3,:)*q_dot+0.5*q_dot.'*Jl2(4:6,:).'*T02(1:3,1:3)*L(2).I*T02(1:3,1:3).'*Jl2(4:6,:)*q_dot;
Tk(3)=0.5*L(3).m*q_dot.'*Jl3(1:3,:).'*Jl3(1:3,:)*q_dot+0.5*q_dot.'*Jl3(4:6,:).'*T03(1:3,1:3)*L(3).I*T03(1:3,1:3).'*Jl3(4:6,:)*q_dot;
Tk(4)=0.5*L(4).m*q_dot.'*Jl4(1:3,:).'*Jl4(1:3,:)*q_dot+0.5*q_dot.'*Jl4(4:6,:).'*T04(1:3,1:3)*L(4).I*T04(1:3,1:3).'*Jl4(4:6,:)*q_dot;
Tk(5)=0.5*L(5).m*q_dot.'*Jl5(1:3,:).'*Jl5(1:3,:)*q_dot+0.5*q_dot.'*Jl5(4:6,:).'*T05(1:3,1:3)*L(5).I*T05(1:3,1:3).'*Jl5(4:6,:)*q_dot;
Tk(6)=0.5*L(6).m*q_dot.'*Jl6(1:3,:).'*Jl6(1:3,:)*q_dot+0.5*q_dot.'*Jl6(4:6,:).'*T06(1:3,1:3)*L(6).I*T06(1:3,1:3).'*Jl6(4:6,:)*q_dot;
Tt=Tk(1)+Tk(2)+Tk(3)+Tk(4)+Tk(5)+Tk(6);



Ul=sym(zeros(6,1));
Ul(1)=-L(1).m*g*Pc1(3);
Ul(2)=-L(2).m*g*Pc2(3);
Ul(3)=-L(3).m*g*Pc3(3);
Ul(4)=-L(4).m*g*Pc4(3);
Ul(5)=-L(5).m*g*Pc5(3);
Ul(6)=-L(6).m*g*Pc6(3);
Ut=Ul(1)+Ul(2)+Ul(3)+Ul(4)+Ul(5)+Ul(6);

La=sym(0);
La=Tt-Ut;


%拉格朗日函数求导
La_q_dot=sym(zeros(6,1));
La_q_dot(1)=diff(La,dt1);
La_q_dot(2)=diff(La,dt2);
La_q_dot(3)=diff(La,dt3);
La_q_dot(4)=diff(La,dt4);
La_q_dot(5)=diff(La,dt5);
La_q_dot(6)=diff(La,dt6);

M=[diff(La_q_dot(1),dt1) diff(La_q_dot(1),dt2) diff(La_q_dot(1),dt3) diff(La_q_dot(1),dt4) diff(La_q_dot(1),dt5) diff(La_q_dot(1),dt6);
   diff(La_q_dot(2),dt1) diff(La_q_dot(2),dt2) diff(La_q_dot(2),dt3) diff(La_q_dot(2),dt4) diff(La_q_dot(2),dt5) diff(La_q_dot(2),dt6);
   diff(La_q_dot(3),dt1) diff(La_q_dot(3),dt2) diff(La_q_dot(3),dt3) diff(La_q_dot(3),dt4) diff(La_q_dot(3),dt5) diff(La_q_dot(3),dt6);
   diff(La_q_dot(4),dt1) diff(La_q_dot(4),dt2) diff(La_q_dot(4),dt3) diff(La_q_dot(4),dt4) diff(La_q_dot(4),dt5) diff(La_q_dot(4),dt6);
   diff(La_q_dot(5),dt1) diff(La_q_dot(5),dt2) diff(La_q_dot(5),dt3) diff(La_q_dot(5),dt4) diff(La_q_dot(5),dt5) diff(La_q_dot(5),dt6);
   diff(La_q_dot(6),dt1) diff(La_q_dot(6),dt2) diff(La_q_dot(6),dt3) diff(La_q_dot(6),dt4) diff(La_q_dot(6),dt5) diff(La_q_dot(6),dt6)];


% C1=[diff(La_q_dot(1),t1) diff(La_q_dot(1),t2) diff(La_q_dot(1),t3) diff(La_q_dot(1),t4) diff(La_q_dot(1),t5) diff(La_q_dot(1),t6);
%     diff(La_q_dot(2),t1) diff(La_q_dot(2),t2) diff(La_q_dot(2),t3) diff(La_q_dot(2),t4) diff(La_q_dot(2),t5) diff(La_q_dot(2),t6);
%     diff(La_q_dot(3),t1) diff(La_q_dot(3),t2) diff(La_q_dot(3),t3) diff(La_q_dot(3),t4) diff(La_q_dot(3),t5) diff(La_q_dot(3),t6);
%     diff(La_q_dot(4),t1) diff(La_q_dot(4),t2) diff(La_q_dot(4),t3) diff(La_q_dot(4),t4) diff(La_q_dot(4),t5) diff(La_q_dot(4),t6);
%     diff(La_q_dot(5),t1) diff(La_q_dot(5),t2) diff(La_q_dot(5),t3) diff(La_q_dot(5),t4) diff(La_q_dot(5),t5) diff(La_q_dot(5),t6);
%     diff(La_q_dot(6),t1) diff(La_q_dot(6),t2) diff(La_q_dot(6),t3) diff(La_q_dot(6),t4) diff(La_q_dot(6),t5) diff(La_q_dot(6),t6)];
% 
% C2=-[diff(La,t1);
%     diff(La,t2);
%     diff(La,t3);
%     diff(La,t4);
%     diff(La,t5);
%     diff(La,t6)];


C3=sym(zeros(6,6));
for e=1:6
    for o=1:6
        C3(e,o)=[diff(M(e,o),q(1)) diff(M(e,o),q(2)) diff(M(e,o),q(3)) diff(M(e,o),q(4)) diff(M(e,o),q(5)) diff(M(e,o),q(6))]*q_dot;
    end
end



%正动力学计算：M(q)q''+C1(q,q')q'+C2(q,q')=tao

%动力学微分方程系数矩阵计算：

G=sym(zeros(6,1));
G=[diff(Ut,t1);
   diff(Ut,t2);
   diff(Ut,t3);
   diff(Ut,t4);
   diff(Ut,t5);
   diff(Ut,t6)];

temp=...
0.5*L(1).m*q_dot.'*Jl1(1:3,:).'*Jl1(1:3,:)+0.5*q_dot.'*Jl1(4:6,:).'*T01(1:3,1:3)*L(1).I*T01(1:3,1:3).'*Jl1(4:6,:)+...
0.5*L(2).m*q_dot.'*Jl2(1:3,:).'*Jl2(1:3,:)+0.5*q_dot.'*Jl2(4:6,:).'*T02(1:3,1:3)*L(2).I*T02(1:3,1:3).'*Jl2(4:6,:)+...
0.5*L(3).m*q_dot.'*Jl3(1:3,:).'*Jl3(1:3,:)+0.5*q_dot.'*Jl3(4:6,:).'*T03(1:3,1:3)*L(3).I*T03(1:3,1:3).'*Jl3(4:6,:)+...
0.5*L(4).m*q_dot.'*Jl4(1:3,:).'*Jl4(1:3,:)+0.5*q_dot.'*Jl4(4:6,:).'*T04(1:3,1:3)*L(4).I*T04(1:3,1:3).'*Jl4(4:6,:)+...
0.5*L(5).m*q_dot.'*Jl5(1:3,:).'*Jl5(1:3,:)+0.5*q_dot.'*Jl5(4:6,:).'*T05(1:3,1:3)*L(5).I*T05(1:3,1:3).'*Jl5(4:6,:)+...
0.5*L(6).m*q_dot.'*Jl6(1:3,:).'*Jl6(1:3,:)+0.5*q_dot.'*Jl6(4:6,:).'*T06(1:3,1:3)*L(6).I*T06(1:3,1:3).'*Jl6(4:6,:);


c=[diff(temp,t1);
   diff(temp,t2);
   diff(temp,t3);
   diff(temp,t4);
   diff(temp,t5);
   diff(temp,t6)];

C=C3-c;

ftip = [1;2;3];  %都是在末端系(6系)的表示  
tautip = [4;5;6];
T06_=double(subs(T06,{t1;t2;t3;t4;t5;t6;dt1;dt2;dt3;dt4;dt5;dt6},{0.1;0.2;0.3;0.4;0.5;0.6;0.6;0.5;0.4;0.3;0.2;0.1}));
ftip_0 = T06_(1:3,1:3) * ftip;
tautip_0 = T06_(1:3,1:3) * tautip;
%动力学微分方程：Mq''+ Cq'+ G + J' * [ftip_0;tautip_0] = tao
J=sym(zeros(6,6));
for i=1:6
    J(:,i)=[cross(z(:,i),P(:,6)-P(:,i));z(:,i)];  %这里的ftip和tautip加在末端系(即6系)的原点
end


m_=double(subs(M,{t1;t2;t3;t4;t5;t6},{0.1;0.2;0.3;0.4;0.5;0.6}));
c_=double(subs(C,{t1;t2;t3;t4;t5;t6;dt1;dt2;dt3;dt4;dt5;dt6},{0.1;0.2;0.3;0.4;0.5;0.6;0.6;0.5;0.4;0.3;0.2;0.1}));


g_=-double(subs(G,{t1;t2;t3;t4;t5;t6;dt1;dt2;dt3;dt4;dt5;dt6},{0.1;0.2;0.3;0.4;0.5;0.6;0.6;0.5;0.4;0.3;0.2;0.1}));
J_ = double(subs(J,{t1;t2;t3;t4;t5;t6;dt1;dt2;dt3;dt4;dt5;dt6},{0.1;0.2;0.3;0.4;0.5;0.6;0.6;0.5;0.4;0.3;0.2;0.1}));


% 正确的tao是这个
% tao=[8.1858;-6.7622;-5.4918;3.4702;7.0105;6.0004]; 




q = [0.1 0.2 0.3 0.4 0.5 0.6];
dq = [0.6 0.5 0.4 0.3 0.2 0.1];
ddq = [4.1 1.2 2.3 3.4 6.1 1.2];
disp('质量矩阵，拉格朗日：');
disp(m_);
disp('质量矩阵，工具箱函数：');
disp(robot.inertia(q));

disp('科氏力与向心力项，拉格朗日：');
disp(c_*dq.');
disp('科氏力与向心力项，工具箱函数：');
disp(robot.coriolis(q, dq)*dq.');

disp('重力项，拉格朗日：');
disp(g_');
disp('重力项，工具箱函数：');
disp(robot.gravload(q));

disp('负载项');
disp((J_.' * [ftip_0;tautip_0])');

tao=m_*ddq.'+c_*dq.'+g_+J_.' * [ftip_0;tautip_0]


inv(m_)*(tao-(c_*dq.'+g_+J_.' * [ftip_0;tautip_0]))

