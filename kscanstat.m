function [Lambda,Start,Stop,Cz,Nz]=kscanstat(x,stepsize)

% calculate Kulldorff's scan statistic
N=length(x); C=length(find(x));
if (C<2), error('C<2!'); end 
I=find(x); LogL=[]; Start=[]; Stop=[]; Cz=[]; Nz=[];
for j=1:length(I)
    % fprintf('\n\nprocessing center position = %d\n',I(j));
    % create zones of increasing "radius" with I(j) as center
    r=0; % start with zero radius, then increment by stepsize
    while 1
        % take zone with radius r
        % fprintf('... r=%d\t',r);
        start=I(j)-r; stop=I(j)+r;
        % fprintf('start=%d,stop=%d\n',start,stop);
        if ((start<1)||(stop>N)), break, end
        z=x(start:stop); nz=length(z);
        if ((nz/N)>0.5), break, end
        zbar=[];
        if ((start-1)>=1), zbar=[zbar, x(1:start-1)]; end % left side of z
        if ((stop+1)<=N), zbar=[zbar, x(stop+1:N)]; end % right side of z
        cz=length(find(z)); p=cz/nz; q=(C-cz)/(N-nz);
        if p==0, p=eps; elseif p==1, p=1-eps; end
        if q==0, q=eps; elseif q==1, q=1-eps; end
        if (p>q)
            logL=cz*log(p)+(nz-cz)*log(1-p)+(C-cz)*log(q)+((N-nz)-(C-cz))*log(1-q);            
            if isnan(logL), error('logL is nan'), end
        else
            logL=C*log(C)+(N-C)*log(N-C)-N*log(N);
            if isnan(logL), error('logL is nan'), end
        end
        LogL=[LogL; logL]; Start=[Start; start]; Stop=[Stop; stop]; Cz=[Cz; cz]; Nz=[Nz; nz];
        r=r+stepsize;
    end    
end
logL0=C*log(C)+(N-C)*log(N-C)-N*log(N);
Lambda=LogL-logL0;
% sort Lambda in descending order, apply to Start and Stop
[Lambda,order]=sort(Lambda,'descend'); Start=Start(order); Stop=Stop(order); Cz=Cz(order); Nz=Nz(order);
