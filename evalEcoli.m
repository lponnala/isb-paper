clear; clc;

load EcoliGBK.mat;
K={'CTA', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACA', 'AGG', 'TTA', 'GTC'}';
selectc=K;
minnumcodons=200; maxnumcodons=800;
L=[]; for i=1:length(S), L=[L; S(i).stop-S(i).start+1]; end
numgenesinthisrange=length(find((L>=(3*minnumcodons))&(L<=(3*maxnumcodons))));
stepsize=1; Nr=100; pvalcrit=0.05;
Laa=[]; Lrc=[]; X={}; H=[]; Nc=[]; Name={}; SignifCz={}; SignifNz={}; SignifD={}; P={}; SS={};
count=0;
for i=1:length(S)
    seq=S(i).sequence;
    if rem(length(seq),3)~=0, continue, end
    numcodons=length(seq)/3;
    if ((numcodons<minnumcodons)||(numcodons>maxnumcodons)), continue, end
    count=count+1;
    fprintf('\nProcessing gene %d of %d...\n',count,numgenesinthisrange);
    if (rem(count,50)==0), fprintf('So far: Number of genes examined = %d, Number of genes with clusters = %d\n',length(H),length(find(H))), end
    x=nan(1,numcodons-1);
    for j=1:numcodons-1 % ignore the last (stop) codon
        codon=seq(3*(j-1)+1:3*j);
        mm=find(strcmpi(codon,selectc));
        if isempty(mm)
            x(j)=0;
        else
            if length(mm)==1
                x(j)=1;
            else
                error('multiple matches!');
            end
        end
    end
    if ~isempty(find(isnan(x))), error('NaNs remain in the x vector!'); end
    if (length(find(x))<2), continue, end
    N=length(x); C=length(find(x));
    Laa=[Laa; N]; Lrc=[Lrc; C];
    fprintf('Number of codons = %d, Number of rare ones = %d\n',N,C);
    X=[X; x];
    tic
    [Lambda,Start,Stop,Cz,Nz]=kscanstat(x,stepsize);
    LambdaMC=kscanstatmc(Nr,N,C,stepsize);
    k=1; pvalest=length(find(LambdaMC>Lambda(k)))/Nr;
    if (pvalest>pvalcrit)
        % fprintf('No clustering!\n');
        H=[H; 0]; Nc=[Nc; 0]; Name=[Name; S(i).name]; SignifCz=[SignifCz; NaN]; SignifNz=[SignifNz; NaN]; SignifD=[SignifD; NaN];
        P=[P; NaN]; SS=[SS; NaN];
    else
        [sLambda,sStart,sStop,sCz,sNz,sP]=mergeclust(Lambda,Start,Stop,Cz,Nz,LambdaMC,Nr,pvalcrit);
        ss=[sStart'; sStop']; sD=sCz./sNz;
        H=[H; 1]; Nc=[Nc; length(sLambda)]; Name=[Name; S(i).name]; P=[P; sP']; SS=[SS; ss];
        SignifCz=[SignifCz; sCz']; SignifNz=[SignifNz; sNz']; SignifD=[SignifD; sD'];
    end
    toc
    save B_Laa.mat Laa; % Laa excludes the last (stop) codon
    save B_Lrc.mat Lrc; % Lrc is the number of rare codons (out of total Laa codons in each gene)
    save B_X.mat X;
    save B_Name.mat Name;
    save B_H.mat H;
    save B_Nc.mat Nc;
    save B_SignifCz.mat SignifCz;
    save B_SignifNz.mat SignifNz;
    save B_SignifD.mat SignifD;
    save B_P.mat P;
    save B_SS.mat SS;
end
% -------------------------------
