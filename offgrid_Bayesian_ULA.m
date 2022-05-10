function H=offgrid_Bayesian_ULA(Y,X,N,N_point,etc)

% N: number of antennas
% N_point: number of grid points
% etc:  number of active grid points

% Note that the ccurrently used grid  uniformly covers the range [-1,1], which is different from the definition ([-90,90]) in Section III-A of TSP2018_Dai . 
% The two forms are equvalent. But the first one can connect with the DFT basis directly (with N=N_point).


[T,M]=size(Y);
search_area=[-1:2/N_point:1];          % the whole grid
reslu=search_area(3)-search_area(2);   % grid interval
F=exp(-1i*pi*(0:N-1)'*search_area)/sqrt(N);
a=0.0001;b=0.0001; 
maxiter=400;
tol=1e-4;

% initialization
converged = false;
iter = 0;
Aw=X*F;
alpha=1;
delta_inv=mean(abs(Aw'*Y), 2)/(norm(Aw))^2;
delta_last=100;

% algrithm begins
while ~converged
       %calculate mu and Sigma
       Phi=Aw;
       Phi_delta = Phi *  diag(delta_inv);
       V_temp= 1/alpha*eye(T) + Phi_delta * Phi';
       Sigma = diag(delta_inv) -Phi_delta' * (V_temp \Phi_delta);
       mu = alpha * Sigma * Phi' * Y;
       switch iter-floor(iter/3)*3
               case 0
                       %update alpha
                       gamma1 = 1 - real(diag(Sigma)) ./ (delta_inv); 
                       resid=Y-Phi*mu;
                       alpha=( M*T + a )/( b +  norm(resid, 'fro')^2  +  M / alpha * sum(gamma1)     );
               case 1
                       %update delta
                       delta_last = delta_inv;

                       sum_mu=sum( mu.*conj(mu), 2);
                       temp=sum_mu + M*real(diag(Sigma));    
                       delta_inv=    ( b+ real(temp)  )/(a+1);
               case 2
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%% grid refine
                        sum_mu=sum( mu.*conj(mu), 2);
                        resid=Y-Phi*mu;
                        Pm=sum_mu;
                        [~,sort_ind]=sort(Pm, 'descend');    
                        index_amp = sort_ind(1:etc);
                        tempPS=Phi*  Sigma(:,index_amp) ;   
                        df=zeros(length(index_amp),1);
                        for j=1:length(index_amp)
                                ii=index_amp(j);
                                ai= exp(-1i*pi*(0:N-1)'* search_area(ii))/sqrt(N);
                                mut=mu(ii,:);
                                Sigmat=Sigma(:,ii);
                                c1=mut*mut' +  M* Sigmat(ii);
                                c1=abs(c1)*(-alpha);
                                Yti=resid +  Phi(:, ii)*mu(ii,:);
                                c2=  M*(  tempPS(:,j) - Phi(:,ii)*Sigmat(ii) )  -Yti*(mut');
                                c2= c2*(-alpha);
                                phii=Phi(:,ii);
                                sinta= search_area(ii); costa=cos(asin(sinta));
                                c3=(-1i*pi* costa)*[0:N-1]';    %  c3=(-1i*2*pi/sqrt(N))*[0:N-1]';
                                tt1=  X*(c3.*ai); 
                                f1= tt1'*phii*c1  +   tt1'*c2;
                                f1= 2*real(f1);
                                df(j)=f1;
                         end    

                         ddff=sign(df.')*reslu/100;                   
                         search_area(index_amp) = search_area(index_amp) +  ddff;
                         F_active=exp(-1i*pi*(0:N-1)'* search_area(index_amp))/sqrt(N);
                         F(:,index_amp)=F_active;
                         Aw(:,index_amp)=X*F_active;
        end


        % stopping criteria
          erro=norm(delta_inv - delta_last)/norm(delta_last);
        if erro < tol || iter >= maxiter
            converged = true;
        end
         iter = iter + 1;   % The definition of "iter" is different from the paper.  Its value should be divided by 3.
   
end


Pm=mean(mu.*conj(mu),2);
[~,sort_ind]=sort(Pm, 'descend');    
ther1=mean(Pm(sort_ind(1:etc)))*0.01;   
index_amp=find(Pm>ther1);
Tn= X*F(:,index_amp);
S= Tn \ Y;
H= F(:,index_amp)*S;










