function[mu_vec,std_vec,Skewness_vec,trend_vec,Times_vec]=give_vecs(CV,Y,cp)
lengthh=length(cp)-1;
T=CV;
mu_vec=zeros(1,lengthh);
std_vec=zeros(1,lengthh);
Skewness_vec=zeros(1,lengthh);
trend_vec=zeros(1,lengthh);
Times_vec=zeros(1,lengthh);
for ii=1:1:lengthh
mu_vec(ii)=mean(T(cp(ii):1:cp(ii+1)));
std_vec(ii)=std(T(cp(ii):1:cp(ii+1)));
Mediann=median(T(cp(ii):1:cp(ii+1)));
%Skewness appriximated by Pearson
%see: https://mathworld.wolfram.com/PearsonModeSkewness.html
Skewness_vec(ii)=(mu_vec(ii)-Mediann)/3*(std_vec(ii));
abc=polyfit(cp(ii):1:cp(ii+1),Y(cp(ii):1:cp(ii+1)),1);
trend_vec(ii)=abc(1);
Times_vec(ii)=cp(ii+1)-cp(ii);
end


       
