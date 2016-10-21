error_em = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
fin_value =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

m = [100,200,300,400,500,1000,2000,3000,4000,5000,10000,20000,30000,40000,50000,100000,200000,300000,400000,500000];

%Declare constants:
miu = 0.05; SK = 1.05; alpha = 2; beta = 0.1; theta = 0.09; x_1 = 1; x_2 = 0.09; true_result = 1.0512;
mc = 1e5; %Monte carlo

for l = 1:20
    n = m(l);
    Y_1_value = zeros(1,mc); Y_2_value = zeros(1,mc); Y_3_value = zeros(1,mc); Payoff_value = zeros(1,mc);
    for p = 1:mc 
        %Step no and T
        T = 1; dt = 1/n; 
        %Random normal arrays:
        Z_1 = randn(1,n); Z_2 = randn(1,n);
        %Brownean arrays:
        B_1 = zeros(1,n); B_2 = zeros(1,n);

        for k = 2:n %Create brownean paths
            B_1(1,k) = B_1(1,k-1) + sqrt(dt) * Z_1(1,k-1);
            B_2(1,k) = B_2(1,k-1) + sqrt(dt) * Z_2(1,k-1);
        end

        %Y1, Y2, Y3 arrays:
        Y_1 = zeros(1,n); Y_2 = zeros(1,n); Y_3 = zeros(1,n); Payoff = zeros(1,n);

        %Calculate Y1, Y2, Y3:
        Y_1(1,1) = x_1; Y_2(1,1) = x_2; Y_3(1,1) = x_1; %ICs

        for k = 2:n
            dw_1 = B_1(1,k) - B_1(1,k-1);
            dw_2 = B_2(1,k) - B_2(1,k-1);
            Y_1(1,k) = Y_1(1,k-1) + miu*Y_1(1,k-1)*dt + Y_1(1,k-1)*sqrt(Y_2(1,k-1))*dw_1;
            Y_2(1,k) = Y_2(1,k-1) + alpha*(theta - Y_2(1,k-1))*dt + beta*sqrt(Y_2(1,k-1))*dw_2;
            Y_3(1,k) = Y_3(1,k-1) + Y_1(1,k-1)*dt;
            Payoff(1,k) = max(Y_3(1,k)/k-SK,0);
        end
        
        Y_1_value(1,p) = Y_1(1,n);
        Y_2_value(1,p) = Y_2(1,n);
        Y_3_value(1,p) = Y_3(1,n);
        Payoff_value(1,p) = Payoff(1,n);
    end
    
    fin_value(l) = mean(Y_1_value);
    error_em(l) = abs(mean(Y_1_value) - true_result);
end


plot(fin_value)



    


