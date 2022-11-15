clear;

L=10;        %number of DoAs
N=100;       % number of antennas
T=60;        % Snapshot
N_point=200; % number of grid points
SNR=10;
Ps=sqrt((   10.^(SNR/10)   )/2);
X=Ps*( randn(T,N)+1i*randn(T,N)  );
DOA=(rand(1,L)-0.5)*2*90;
A = exp(-1i*pi*(0:N-1)'*sind(DOA));
H=A*(randn(L,1)+1i*randn(L,1));
noise=sqrt(1/2)*(randn(T,1)+1i*randn(T,1));
Y=X*H + noise;
etc=T;       % or etc=length(DOA) if the number of DoAs is known

%% Proposed method
H_est=offgrid_Bayesian_ULA(Y,X,N,N_point,etc);
norm(H-H_est,'fro')^2/norm(H,'fro')^2
