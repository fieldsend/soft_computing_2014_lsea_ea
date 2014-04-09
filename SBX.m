function [x_c, x_d] = SBX(x1,x2,SBX_n,mn,mx)

% simulated binary crossover

l = length(x1);
x_c =x1;
x_d = x2;
for k=1:l
    if rand()<0.5 %0.5 probability of crossing over variable
        valid=0;
        while (valid==0)
            u=rand();
            if (u<0.5)
                beta=(2*u)^(1/(SBX_n+1));
            else
                beta=(0.5/(1-u))^(1/(SBX_n+1));
            end
            mean_p=0.5*(x1(k)+x2(k));
            c1=mean_p-0.5*beta*abs(x1(k)-x2(k));
            c2=mean_p+0.5*beta*abs(x1(k)-x2(k));
            if (c1>=mn(k) && c1<=mx(k)) && (c2>=mn(k) && c2<=mx(k))
                valid=1;
            end
        end
        x_c(k) = c1;
        x_d(k) = c2;
    end
end