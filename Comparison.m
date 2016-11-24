clc
clear
close all

%% Data I/O %%
dir1 = 'DATA'
dir2 = 'RESULT'

fn1 = 'CHANNEL_0180_profile.plt'
fn2 = 'U.plt'
fn3 = 'k.plt'
fn4 = 'dissipation.plt'

strcat(dir1,'/',fn1)

data_p1 = fopen(strcat(dir1,'/',fn1));
data_p2 = fopen(strcat(dir2,'/',fn2));
data_p3 = fopen(strcat(dir2,'/',fn3));
data_p4 = fopen(strcat(dir2,'/',fn4));

A = fscanf(data_p1,'%f %f %f %f %f %f %f %f',[8 inf]);
B = fscanf(data_p2,'%f %f %f %f %f',[5 inf]);
C = fscanf(data_p3,'%f %f %f %f',[4 inf]);
D = fscanf(data_p4,'%f %f %f %f',[4 inf]);

fclose(data_p1)
fclose(data_p2)
fclose(data_p3)
fclose(data_p4)

%% Substitution %%
Y   = B(1,:);
Yp  = B(2,:);
U   = B(4,:);
k   = C(4,:);
dis = D(4,:);

Ny = size(B);
Ny = Ny(2);

Y_exac   = A(1,:);
Yp_exac  = A(2,:);
U_exac   = A(4,:);
k_exac   = A(6,:);
dis_exac = A(8,:);

%% Plotting %%
figure
semilogx(Yp(1:Ny/2),U(1:Ny/2))
hold on
semilogx(Yp_exac,U_exac,'o')

figure
plot(Yp(1:Ny/2),k(1:Ny/2))
hold on
plot(Yp_exac,k_exac,'o')

figure
plot(Yp(1:Ny/2),dis(1:Ny/2))
hold on
plot(Yp_exac,dis_exac,'o')
