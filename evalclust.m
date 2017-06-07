
% //// This script answers the questions ////
% NOTE: the word "cluster" means a statistically significant cluster of
% rare codons, found by the ScanStat method of Kulldorff (1997)

% -----------------------------------------------
% mRNA-protein Vs cluster-density
% -----------------------------------------------
load ChurchData.mat;
% columns contain, in order:
% name,start,stop,strand,func,seq,signal,nabd,spint,rnacopies/cell
name=D(:,1); start=cell2mat(D(:,2)); stop=cell2mat(D(:,3)); strand=D(:,4); func=D(:,5);
seq=D(:,6); signal=D(:,7);
prot=cell2mat(D(:,8)); % protein (N-abd in molecules/cell)
spint=cell2mat(D(:,9)); % rna "SP Intensity"
rna=cell2mat(D(:,10)); % rna (estimated copies/cell)
ratio=prot./rna; hist(ratio), title('histogram of prot/rna ratio')

load B_Name.mat; load B_H.mat; load B_Nc.mat;
load B_SignifNz.mat; load B_SignifCz.mat; load B_SignifD.mat;
myname={}; myratio=[]; myprot=[]; myrna=[];
myNc=[]; myNz=[]; myCz=[]; myD=[];
for i=1:length(name)
    ind=strmatch(name{i},Name,'exact');
    if isempty(ind), continue, end
    if length(ind)>1, error('multiple matches for %s\n',name{i}); end    
    if Nc(ind)==0
        % myNc=[myNc; 0];
        % myNz=[myNz; NaN]; myD=[myD; NaN];
    else
        myname=[myname; name{i}]; myratio=[myratio; ratio(i)]; myprot=[myprot; prot(i)]; myrna=[myrna; rna(i)];
        myNc=[myNc; Nc(ind)];
        myNz=[myNz; SignifNz{ind}]; myCz=[myCz; SignifCz{ind}]; myD=[myD; SignifD{ind}];
    end
end

% corr between ratio and cluster-density
[R_D,P_D]=corrcoef(myratio,myD), fprintf('Corr between ratio, D = %f, p-value = %f\n',R_D(1,2),P_D(1,2)); 
figure, plot(myratio,myD,'b*'), xlabel('protein yield'), ylabel('cluster density')

% partial corr between prot and cluster-density, controlling for rna
[p_rho,p_pval]=partialcorr(myprot,myD,myrna)
fprintf('partial corr between prot and cluster-density (while controlling for rna) = %f, p-value = %f\n',p_rho,p_pval)
% -----------------------------------------------

% ---------------------------------------------
% proportion of genes that contain rare codons
% ---------------------------------------------
clear;
load B_Laa.mat; load B_Lrc.mat;
fprintf('percentage of genes that contain atleast one rare codon = %f\n',100*length(find(Lrc))/length(Lrc));
rcprop=100*Lrc./Laa;
figure, hist(rcprop), title('rare codon percentage per gene')
fprintf('on an average, %f percent of each gene is made up of rare codons\n\n\n',mean(rcprop));

% ------------------------------------------------------------------
% proportion of genes that contain doublets/triplets of rare codons
% ------------------------------------------------------------------
clear;
load B_X.mat;
D=[];
for i=1:length(X)
    x=X{i}; I=find(x);
    for j=1:length(I)
        if I(j)~=length(x)
            if x(I(j)+1)==1, D=[D; i]; break, end
        end
    end
end
fprintf('percentage of genes that contain atleast one doublet = %f\n',100*length(D)/length(X));
F=[];
for i=1:length(X)
    x=X{i}; I=find(x);
    for j=1:length(I)
        if ((I(j)>=2)&&(I(j)<=(length(x)-1)))
            xp=x(I(j)-1:I(j)+1);
            if length(find(xp))>=2, F=[F; i]; break, end
        end
    end
end
fprintf('percentage of genes that contain atleast one 2-or-more/3 region = %f\n',100*length(F)/length(X));
T=[];
for i=1:length(X)
    x=X{i}; I=find(x);
    for j=1:length(I)
        if ((I(j)~=1)&&(I(j)~=length(x)))
            if ((x(I(j)-1)==1)&&(x(I(j)+1)==1)), T=[T; i]; break, end
        end
    end
end
fprintf('percentage of genes that contain atleast one triplet = %f\n',100*length(T)/length(X));
F4=[]; F3=[];
for i=1:length(X)
    x=X{i}; I=find(x);
    for j=1:length(I)
        if ((I(j)>=3)&&(I(j)<=(length(x)-2)))
            xp=x(I(j)-2:I(j)+2);
            if length(find(xp))>=4, F4=[F4; i]; break, end
            if length(find(xp))>=3, F3=[F3; i]; break, end
        end
    end
end
fprintf('percentage of genes that contain atleast one 4-or-more/5 region = %f\n',100*length(F4)/length(X));
fprintf('percentage of genes that contain atleast one 3-or-more/5 region = %f\n',100*length(F3)/length(X));
F=[];
for i=1:length(X)
    x=X{i}; I=find(x);
    for j=1:length(I)
        if ((I(j)>=4)&&(I(j)<=(length(x)-3)))
            xp=x(I(j)-3:I(j)+3);
            if length(find(xp))>=6, F=[F; i]; break, end
        end
    end
end
fprintf('percentage of genes that contain atleast one 6-or-more/7 region = %f\n\n\n',100*length(F)/length(X));

% --------------------------------------------------
% percentage of genes that have atleast one cluster
% --------------------------------------------------
clear;
load B_H.mat;
fprintf('percentage of genes that have atleast one cluster = %f\n\n\n',100*length(find(H==1))/length(H));

% --------------------------------------------------
% some cluster-size evaluation metrics
% --------------------------------------------------
clear;
load B_H.mat; load B_Nc.mat; load B_Laa.mat; load B_SignifCz.mat; load B_SignifNz.mat; load B_SignifD.mat; load B_Name.mat;
nz=[]; cz=[]; d=[]; gene={};
for i=1:length(H)
    if (H(i)==1)
        signifNz=SignifNz{i}; signifCz=SignifCz{i}; signifD=SignifD{i};
        for j=1:Nc(i)
            nz=[nz; signifNz(j)];
            cz=[cz; signifCz(j)];
            d=[d; signifD(j)];
            gene=[gene; Name{i}];
        end
    end
end
% max cluster size
fprintf('max cluster size = %d\n',max(nz));
N=find(nz==max(nz)); fprintf('number of clusters having max size = %d\n',length(N));
fprintf('gene(s) having maximum-size clusters:\n'); disp(gene(N))
% max cluster content (i.e. number of rare codons)
fprintf('max cluster content = %d\n',max(cz));
C=find(cz==max(cz)); fprintf('number of clusters having max content = %d\n',length(C));
fprintf('gene(s) having maximum-content clusters:\n'); disp(gene(C))
% max cluster density
fprintf('max cluster density = %d\n',max(d));
D=find(d==max(d)); fprintf('number of clusters having max density = %d\n',length(D));
% fprintf('gene(s) having maximum-density clusters:\n'); disp(gene(D))
fprintf('\n\n');

% --------------------------------------------------
% gap between clusters
% --------------------------------------------------
clear;
load B_Nc.mat; load B_SS.mat; load B_Name.mat;
gap=[]; gene={};
for i=1:length(Nc)
    if (Nc(i)>1)
        ss=SS{i}; start=ss(1,:); stop=ss(2,:);
        [start,I]=sort(start); stop=stop(I);
        for j=2:length(start)
            if (start(j)>stop(j-1))
                gap=[gap; (start(j)-stop(j-1)-1)];
                gene=[gene; Name{i}];
            end
        end
    end
end
figure, hist(gap), title('gap between non-overlapping clusters')
fprintf('Average gap (in codons) between non-overlapping clusters = %f\n\n\n',mean(gap));

% --------------------------------------------------
% percentage of all clusters in short proteins
% --------------------------------------------------
clear;
load B_Laa.mat; load B_H.mat; load B_Nc.mat;
% lothresh=300; I=find(Laa<=lothresh); fprintf('percentage of short (<=%d aa) proteins that have atleast one cluster = %f\n\n\n',lothresh,100*length(find(H(I)))/length(I));
% lothresh=[250 350]; I2=find((Laa>=lothresh(1))&(Laa<=lothresh(2))); fprintf('percentage of short (between %d and %d aa) proteins that have atleast one cluster = %f\n\n\n',lothresh(1),lothresh(2),100*length(find(H(I2)))/length(I2));
% hithresh=500; J=find(Laa>hithresh); fprintf('percentage of long (>%d aa) proteins that have atleast one cluster = %f\n\n\n',hithresh,100*length(find(H(J)))/length(J));
lothresh=300; I=find(Laa<=lothresh); fprintf('percentage of all clusters that occur in short (<=%d aa) proteins = %f\n',lothresh,100*sum(Nc(I))/sum(Nc));
lothresh=[250 350]; I2=find((Laa>=lothresh(1))&(Laa<=lothresh(2))); fprintf('percentage of all clusters that occur in short (between %d and %d aa) proteins = %f\n',lothresh(1),lothresh(2),100*sum(Nc(I2))/sum(Nc));
hithresh=500; J=find(Laa>hithresh); fprintf('percentage of all clusters that occur in long (>%d aa) proteins = %f\n\n\n',hithresh,100*sum(Nc(J))/sum(Nc));

% --------------------------------------------------
% number of proteins having atleast one slow-translating cluster Vs. length
% of protein
% --------------------------------------------------
clear;
load B_Laa.mat; load B_H.mat;
I=find(H==1); L=Laa(I); figure, hist(L), title('length of proteins that contain atleast one cluster')
% disp(mean(L))
thresh=300;
fprintf('Out of %d proteins having atleast one cluster, %d occur among proteins of length>=%d\n\n\n',length(L),length(find(L>=thresh)),thresh);

% --------------------------------------------------
% the positional distribution of slow-translating codons
% --------------------------------------------------
% //// codon analysis ////
clear;
load EcoliGBK.mat;
selectc={'CTA', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACA', 'AGG', 'TTA', 'GTC'}';
minnumcodons=200; maxnumcodons=800;
L=[]; for i=1:length(S), L=[L; S(i).stop-S(i).start+1]; end
numgenesinthisrange=length(find((L>=(3*minnumcodons))&(L<=(3*maxnumcodons))));
X={}; count=0;
for i=1:length(S)
    seq=S(i).sequence;
    if rem(length(seq),3)~=0, continue, end
    numcodons=length(seq)/3;
    if ((numcodons<minnumcodons)||(numcodons>maxnumcodons)), continue, end
    count=count+1;
    % fprintf('Processing %d of %d...\n',count,numverifgenesinthisrange);
    % fprintf('Processing %d of %d...\n',count,numgenesinthisrange);
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
    X{i}=x;
end
% check how many slow-translating codons occur in the first 50 codons
N=0; N50=0;
for i=1:length(X)
    I=find(X{i});
    N=N+length(I);
    N50=N50+length(find(I<=50));
end
fprintf('%f percent of all rare codons occur near translation-starts (in the first 50 codons)\n',100*(N50/N));

% //// cluster analysis ////
clear;
load B_H.mat; load B_SS.mat; load B_Name.mat;
S=[]; gene={};
for i=1:length(H)
    if (H(i)==1)
        ss=SS{i};
        for j=1:size(ss,2)
            S=[S; [ss(1,j),ss(2,j)]];
            gene=[gene; Name{i}];
        end
    end
end
% what percentage of the clusters are completely inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that are completely inside the first 50 amino acids = %f\n',100*length(find(S(:,2)<=50))/size(S,1));
% what percentage of the clusters have atleast one-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast one-aa inside the first 50 amino acids = %f\n',100*length(find(S(:,1)<=50))/size(S,1));
% what percentage of the clusters have atleast 5-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast 5-aa inside the first 50 amino acids = %f\n',100*length(find((S(:,1)<=46)&(S(:,2)>=50)))/size(S,1));
% what percentage of the clusters have atleast 10-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast 10-aa inside the first 50 amino acids = %f\n',100*length(find((S(:,1)<=41)&(S(:,2)>=50)))/size(S,1));
% how many of the clusters are atleast half-inside the first 50 amino acids (150 bases)?
count=0; fcount=0;
for i=1:size(S,1)
    if S(i,2)<=50, fcount=fcount+1; end
    midpointindex=S(i,1)+floor((S(i,2)-S(i,1))/2);
    if (midpointindex<=50), count=count+1; end
end
fprintf('percentage of clusters that are atleast half-inside the first 50 amino acids = %f\n',100*count/size(S,1));
fprintf('percentage of clusters that are completely inside the first 50 amino acids = %f\n',100*fcount/size(S,1));
