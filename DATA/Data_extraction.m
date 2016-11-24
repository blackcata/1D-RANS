clc
clear
close all 

filename = 'CHANNEL_0180_'

data1 = 'mean_prof'
data2 = 'vel_fluc_prof'
data3 = 'balance_budget'

fn1  = strcat(filename,data1,'.plt')
fn2  = strcat(filename,data2,'.plt')
fn3  = strcat(filename,data3,'.plt')

data_p1 = fopen(fn1);
data_p2 = fopen(fn2);
data_p3 = fopen(fn3);

A = fscanf(data_p1,'%f %f %f %f %f %f',[6 inf]); 
B = fscanf(data_p2,'%f %f %f %f %f %f %f %f %f',[9 inf]); 
C = fscanf(data_p3,'%f %f %f %f %f %f %f %f %f',[9 inf]); 

fclose(data_p1)
fclose(data_p2)
fclose(data_p3)

%% Data Transfer %%
Y    = A(1,:);
Yp   = A(2,:);
U    = A(3,:);
DUDY = A(4,:);
W    = A(5,:);
P    = A(6,:);

uu   = B(3,:);
vv   = B(4,:);
ww   = B(5,:);
uv   = B(6,:);
uw   = B(7,:);
vw   = B(8,:);
k    = B(9,:);

prod = C(3,:);
e    = C(8,:);

%% Plotting %%
semilogx(Yp,U,'ro')
xlabel('y+')
ylabel('u+')
axis([1 max(Yp)+20 0 max(U)+2])

%% Constants %%

Re_tau = 180
nu = 3.5e-4
del = 1.0

u_tau = Re_tau * nu /del

k_max = max(k*u_tau*u_tau) 
e_max = max(e)

%% Writing Y,Yp,U,k %%
fileID = fopen('CHANNEL_0180_profile.plt','w');
fprintf(fileID,'VARIABLES="Y","Y+","U","U+","k","k+","diss","diss+" \n');
fprintf(fileID,'%f %f %f %f %f %f %f %f \n',[Y' Yp' U'*u_tau U' k'*u_tau^2 k' e'*u_tau^2 e']');
fclose(fileID);
