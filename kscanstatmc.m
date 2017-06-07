
function LambdaMC=kscanstatmc(Nr,N,C,stepsize)

if (C<2), error('C<2!'); end
LambdaMC=nan(1,Nr);
for k=1:Nr
    v=zeros(1,N);
    ind=randint(1,C,[1 N]);
    v(ind)=1;
    while length(find(v))<C
        % add as many are less
        ind=randint(1,C-length(find(v)),[1 N]);
        v(ind)=1;
    end
    if (length(find(v))~=C), error('incorrect number of points chosen'), end
    LambdaMC(k)=kscanstatmax(v,stepsize);
end
if ~isempty(find(isnan(LambdaMC))), error('some of the LambdaMC values are NaN!'); end
