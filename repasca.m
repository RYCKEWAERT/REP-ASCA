function [model] = rep_asca( X, design, X_rep, d_rep,klimit)
% REP-ASCA
% RYCKEWAERT MAXIME (02/12/2020)
% Reduction of repeatability error for analysis of variance-Simultaneous Component Analysis (REP-ASCA)
% Analysis of variance developed for multivariate dataset
% (included in a design of experiments, experimental design)
%
%example:
% rep_asca( X, design, X_rep, d_rep,klimit)
% X : multivariate dataset 
% design : experimental design (binary : 0/1)
% opts : options for asca and orthogonalisation
% X_rep: repeated measures 
% d_rep: samples identification 
% klimit: limitation on component to analyse
%
%see also :
% M. Ryckewaert, N. Gorretta, F. Henriot, F. Marini, et J.-M. Roger,
% « Reduction of repeatability error for analysis of variance-Simultaneous Component Analysis (REP-ASCA): Application to NIR spectroscopy on coffee sample »,
% Analytica Chimica Acta, vol. 1101, p. 23?31, mars 2020, doi: 10.1016/j.aca.2019.12.024.

if nargin==2
     model = asca_final( X, design);
    return;
end;


BS =  d_rep* (inv(d_rep'*d_rep))*d_rep'*X_rep;
WS = X_rep - BS;
[~,~,L_err]=svds((WS),klimit);

results_asca = asca_min(X,design);
list_factor = results_asca.TermLabels(2:end);
for j=1:length(list_factor)
    eval(sprintf(' explained_var(1,j) = results_asca.X%s.EffectExplVar ;',strrep(list_factor{j},' ','') ));
    list_factor_name{j} = strcat('X',strrep(list_factor{j},' ',''));
end
explained_var(1,size(list_factor,2)+1) = results_asca.XRes.EffectExplVar;
list_factor_name{size(list_factor,2)+1} = 'XRes';

for i = 1:size(L_err,2)
    k_W = L_err(:,1:i)';
    X_bar = X-X*k_W'*k_W;
    results_asca = asca_min(X_bar,design);
    list_factor = results_asca.TermLabels(2:end);
    for j=1:length(list_factor)
        eval(sprintf(' explained_var(i+1,j) = results_asca.X%s.EffectExplVar ;',strrep(list_factor{j},' ','') ));
        list_factor_name{j} = strcat('X',strrep(list_factor{j},' ',''));
    end
    explained_var(i+1,size(list_factor,2)+1) = results_asca.XRes.EffectExplVar;
    list_factor_name{size(list_factor,2)+1} = 'XRes';
end

model.explained_var= explained_var;
model.list_factor_name= list_factor_name;
model.L_err = L_err;

end 

%% ASCA 
function [model] = asca_final( X, design)


dmatr = designmat_fct(design);
[dmatr_factor,terms,labels, order]=design_creator(dmatr);

Xcomputed=X;

%Work by factor: 
for i=1:length(dmatr_factor);
    X_i=dmatr_factor{i}*pinv(dmatr_factor{i})*Xcomputed;
    ssq_i=sum(sum(X_i.^2));
    Xcomputed=Xcomputed-X_i;
    l=strtrim(labels{i});
    eval(['model.X', l,'.EffectMatrix=X_i;'])
    eval(['model.X', l,'.EffectSSQ=ssq_i;'])
    
    if i==1
        ssqtot=sum(sum(Xcomputed.^2));
        model.Xdata.CenteredData=Xcomputed;
        model.Xdata.CenteredSSQ=ssqtot;
    else
        expVar=100*(ssq_i/ssqtot);
        eval(['model.X', l,'.EffectExplVar=expVar;'])
    end
    
    eval(['model.X', l,'.DesignMatrix=dmatr_factor{i};'])
    eval(['model.X', l,'.DesignTerms=terms{i};'])
    eval(['model.X', l,'.EffectLabel=strtrim(labels{i});'])
    eval(['model.X', l,'.TermOrder=order(i);'])
end


for i=2:length(dmatr_factor);
    l=strtrim(labels{i});
   
    X_i=model.Xdata.CenteredData;
    for j=1:length(l);
        eval(['X_i=X_i-model.X',l,'.EffectMatrix;'])
    end
    eval(['model.X', l,'.ReducedMatrix=X_i;'])
end


model.XRes.EffectMatrix=Xcomputed;
model.XRes.EffectSSQ=sum(sum(Xcomputed.^2));
model.XRes.EffectExplVar=100*(model.XRes.EffectSSQ/ssqtot);
model.XRes.EffectLabel='Res';
model.XRes.TermOrder=max(order)+1;
model.TermLabels=labels;

model=permutation_test(model);
%SCA Models
model=scastep(model);

end


%% Permutation test
function newmodel=permutation_test(model)
nperm = 500;
newmodel=model;
labels=model.TermLabels;
signfacts=cell(length(labels)-1,1);
sc=0;
for i=2:length(labels)
    l=strtrim(labels{i}); 
    eval(['Xr=model.X', l, '.ReducedMatrix; '])
    eval(['Dr=model.X', l, '.DesignMatrix; '])
    ssqp=ptest(Xr,Dr, 500);
    eval(['seff=model.X', l, '.EffectSSQ; '])
    p=length(find(ssqp>=seff))./nperm;
    if p<=0.05
        sc=sc+1;
        signfacts{sc}=l;
    end  
    eval(['newmodel.X',l,'.EffectSignif.NullDistr=ssqp;'])
    eval(['newmodel.X',l,'.EffectSignif.p=p;'])
end
signfacts=signfacts(1:sc);
newmodel.SignificantTerms=signfacts;
end
%% Design creator
function [dmatr,terms, labels, order]=design_creator(design)
nfactor=length(design);
nsize=size(design{1},1);

indmat=fullfact(repmat(2,1,nfactor));
nmat=size(indmat,1);
dmatr=cell(1,nmat);
terms=cell(1,nmat);
labels=cell(1,nmat);
order=zeros(nmat,1);

for i=1:nmat
    Dm=ones(nsize,1);
    for j=1:nfactor
        if indmat(i,j)==1
            effmat=design{j};
        else
            effmat=ones(nsize,1);
        end
        Dm=kron(Dm,effmat);
        Dm=Dm(1:nsize+1:end,:);
    end
    dmatr{i}=Dm;
    terms{i}=find(indmat(i,:)==1);
    labels{i}=char(64+find(indmat(i,:)==1));
    order(i)=length(terms{i});   
end
labels=sort(char(labels),2);
%Sorting according to increasing order of interactions
[labels, newindex]=sortrows(labels);

terms=terms(newindex);
dmatr=dmatr(newindex);
order=order(newindex);
labels=cellstr(labels)';
labels{1}='Mean';
end
function ssqp=ptest(X,D, nperm)

ns=size(X,1);  %Number of samples
ssqp=zeros(nperm,1);   %Initialization of the permuted SSQ vector

for i=1:nperm
    hh=randperm(ns);
    Xpp=D(hh,:)*pinv(D(hh,:))*X;
    ssqp(i)=sum(sum(Xpp.^2));
end
end
function dmatr=designmat_fct(design)
nsize=size(design,1);
nfactor=size(design,2);
dmatr=cell(1,nfactor);
for i=1:nfactor
    lev=unique(design(:,i));
    nl=length(lev);
    dmat=zeros(nsize,nl-1);
    for j=1:nl-1
        dmat(design(:,i)==lev(j),j)=1;
    end
    dmat(design(:,i)==lev(nl),:)=-1;
    dmatr{i}=dmat;
end
end

%% SCA
function smodel=scastep(ascamodel)

smodel=ascamodel;
dlab=ascamodel.TermLabels;

for i=2:length(dlab)
    l=strtrim(dlab{i});
    eval(['Xr=ascamodel.X', l, '.EffectMatrix; '])
    eval(['ssqr=ascamodel.X', l, '.EffectSSQ; '])  
    R=rank(Xr);
    [u,s,P]=svds(Xr,R);
    T=u*s;
    explainedVar=100*(diag(s).^2)/ssqr;
    eval(['smodel.X', l, '.SCA.Model.Scores=T; '])
    eval(['smodel.X', l, '.SCA.Model.Loadings=P; '])
    eval(['smodel.X', l, '.SCA.Model.ExplVar=explainedVar; '])
end
end

function [model] = asca_min( X, design)


dmatr = designmat_fct(design);
[dmatr_factor,terms,labels, order]=design_creator(dmatr);

Xcomputed=X;

%Work by factor: 
for i=1:length(dmatr_factor);
    X_i=dmatr_factor{i}*pinv(dmatr_factor{i})*Xcomputed;
    ssq_i=sum(sum(X_i.^2));
    Xcomputed=Xcomputed-X_i;
    l=strtrim(labels{i});
    eval(['model.X', l,'.EffectMatrix=X_i;'])
    eval(['model.X', l,'.EffectSSQ=ssq_i;'])
    
    if i==1
        ssqtot=sum(sum(Xcomputed.^2));
        model.Xdata.CenteredData=Xcomputed;
        model.Xdata.CenteredSSQ=ssqtot;
    else
        expVar=100*(ssq_i/ssqtot);
        eval(['model.X', l,'.EffectExplVar=expVar;'])
    end
    
    eval(['model.X', l,'.DesignMatrix=dmatr_factor{i};'])
    eval(['model.X', l,'.DesignTerms=terms{i};'])
    eval(['model.X', l,'.EffectLabel=strtrim(labels{i});'])
    eval(['model.X', l,'.TermOrder=order(i);'])
end

model.XRes.EffectMatrix=Xcomputed;
model.XRes.EffectSSQ=sum(sum(Xcomputed.^2));
model.XRes.EffectExplVar=100*(model.XRes.EffectSSQ/ssqtot);
model.XRes.EffectLabel='Res';
model.XRes.TermOrder=max(order)+1;
model.TermLabels=labels;
end
