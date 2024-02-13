function distance = EgoDist(A1,A2,delta,cap,DistanceType)

%   Reference: 
%   C. Piccardi, Metrics for network comparison using egonet feature 
%   distribution, Scientific Reports, 13, 14657, 2023. 
%   [https://doi.org/10.1038/s41598-023-40938-4]
%
%   ------------------------------------------------------------------
%   INPUTS:
%     
%   A1,A2: Binary undirected adjacency matrices 
%          (possibly with different size)
%
%   delta: discretization interval for discrete distributions (0<delta<1)
%
%   cap: upper bound for discrete distributions (0<cap<=1, cap>>delta)
%
%   DistanceType: {'D','C','P','SUM','CP','DC','DP','DCP'}
%
%   OUTPUT:     
%
%   distance: distance between networks A1, A2
%   ------------------------------------------------------------------
%
%   Example of usage (the two nets are in the folder "networks" of the 
%   zip distribution file):
%
%   load('net_SFBA_n1000_d000_1.mat'); A1=A; %loading A1
%   load('net_GEO_n2000_d000_10.mat'); A2=A; %loading A2
%   distance=EgoDist(A1,A2,0.01,1,'DCP')
%
%   distance =
%
%     672.5809
%   
% 
%   Note: To speed up computations, no check is performed on the
%   correctness and consistency of the inputs.
%   


N1=length(A1);
N2=length(A2);
r=round(cap/delta);        

%preparing egonet statistics

%normalized degree
if strcmp(DistanceType,'D')||strcmp(DistanceType,'SUM')||strcmp(DistanceType,'DC')||strcmp(DistanceType,'DP')||strcmp(DistanceType,'DCP')
    
    deg1=sum(A1);
    dd1=(deg1-min(deg1))/(max(deg1)-min(deg1));
    dd1(isnan(dd1))=0;

    deg2=sum(A2);
    dd2=(deg2-min(deg2))/(max(deg2)-min(deg2));
    dd2(isnan(dd2))=0;

end

%clustering coefficient
if strcmp(DistanceType,'C')||strcmp(DistanceType,'SUM')||strcmp(DistanceType,'DC')||strcmp(DistanceType,'CP')||strcmp(DistanceType,'DCP')
    
    cc1=zeros(1,N1);
    for u=1:N1
        V=find(A1(u,:));
        k=length(V);
        S=A1(V,V);
        cc1(u)=sum(S(:))/(k^2-k);
    end
    cc1(isnan(cc1))=0;

    cc2=zeros(1,N2);
    for u=1:N2
        V=find(A2(u,:));
        k=length(V);
        S=A2(V,V);
        cc2(u)=sum(S(:))/(k^2-k);
    end
    cc2(isnan(cc2))=0;

end

%egonet persistence
if strcmp(DistanceType,'P')||strcmp(DistanceType,'SUM')||strcmp(DistanceType,'CP')||strcmp(DistanceType,'DP')||strcmp(DistanceType,'DCP')
    
    pp1=zeros(1,N1);
    deg1=sum(A1);
    for v=1:N1
        set=logical(A1(v,:)); set(v)=1;
        pp1(v)=sum(sum(A1(set,set)))/sum(deg1(set));
    end     
    pp1(isnan(pp1))=0;

    pp2=zeros(1,N2);
    deg2=sum(A2);
    for v=1:N2
        set=logical(A2(v,:)); set(v)=1;
        pp2(v)=sum(sum(A2(set,set)))/sum(deg2(set));
    end     
    pp2(isnan(pp2))=0;
    
end

%computing distributions and distance

switch DistanceType
    case 'D'
        pdf=zeros(r,1);
        dd1(dd1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(dd1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        D1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        dd2(dd2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(dd2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        D2(:)=cumsum(pdf)/N2;
        
        distance=norm(D1-D2);

    case 'C'
        pdf=zeros(r,1);
        cc1(cc1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(cc1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        C1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        cc2(cc2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(cc2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        C2(:)=cumsum(pdf)/N2;

        distance=norm(C1-C2);

    case 'P'
        pdf=zeros(r,1);
        pp1(pp1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(pp1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        P1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        pp2(pp2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(pp2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        P2(:)=cumsum(pdf)/N2;

        distance=norm(P1-P2);

    case 'SUM'
        pdf=zeros(r,1);
        dd1(dd1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(dd1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        D1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        dd2(dd2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(dd2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        D2(:)=cumsum(pdf)/N2;

        pdf=zeros(r,1);
        cc1(cc1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(cc1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        C1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        cc2(cc2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(cc2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        C2(:)=cumsum(pdf)/N2;

        pdf=zeros(r,1);
        pp1(pp1==1)=1-delta/2; %to put ones in the last bin
        for n=1:N1
            a=floor(pp1(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        P1(:)=cumsum(pdf)/N1;

        pdf=zeros(r,1);
        pp2(pp2==1)=1-delta/2; %to put ones in the last bin
        for n=1:N2
            a=floor(pp2(n)/delta)+1;
            if a<=r
                pdf(a)=pdf(a)+1;
            end
        end
        P2(:)=cumsum(pdf)/N2;

        distance=norm(D1-D2)+norm(C1-C2)+norm(P1-P2);

    case 'DC'
        pdf=zeros(r,r);
        dd1(dd1==1)=1-delta/2; cc1(cc1==1)=1-delta/2;
        for n=1:N1
            a=floor(dd1(n)/delta)+1;
            b=floor(cc1(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        DC1(:,:)=cumsum(cumsum(pdf,1),2)/N1;

        pdf=zeros(r,r);
        dd2(dd2==1)=1-delta/2; cc2(cc2==1)=1-delta/2;
        for n=1:N2
            a=floor(dd2(n)/delta)+1;
            b=floor(cc2(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        DC2(:,:)=cumsum(cumsum(pdf,1),2)/N2;

        distance=norm(DC1-DC2,'fro');
        
    case 'DP'
        pdf=zeros(r,r);
        dd1(dd1==1)=1-delta/2; pp1(pp1==1)=1-delta/2;
        for n=1:N1
            a=floor(dd1(n)/delta)+1;
            b=floor(pp1(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        DP1(:,:)=cumsum(cumsum(pdf,1),2)/N1;

        pdf=zeros(r,r);
        dd2(dd2==1)=1-delta/2; pp2(pp2==1)=1-delta/2;
        for n=1:N2
            a=floor(dd2(n)/delta)+1;
            b=floor(pp2(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        DP2(:,:)=cumsum(cumsum(pdf,1),2)/N2;

        distance=norm(DP1-DP2,'fro');
        
    case 'CP'
        pdf=zeros(r,r);
        cc1(cc1==1)=1-delta/2; pp1(pp1==1)=1-delta/2;
        for n=1:N1
            a=floor(cc1(n)/delta)+1;
            b=floor(pp1(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        CP1(:,:)=cumsum(cumsum(pdf,1),2)/N1;

        pdf=zeros(r,r);
        cc2(cc2==1)=1-delta/2; pp2(pp2==1)=1-delta/2;
        for n=1:N2
            a=floor(cc2(n)/delta)+1;
            b=floor(pp2(n)/delta)+1;
            if a<=r && b<=r
                pdf(a,b)=pdf(a,b)+1;
            end
        end
        CP2(:,:)=cumsum(cumsum(pdf,1),2)/N2;

        distance=norm(CP1-CP2,'fro');

    case 'DCP'
        pdf=zeros(r,r,r);
        cc1(cc1==1)=1-delta/2; pp1(pp1==1)=1-delta/2; dd1(dd1==1)=1-delta/2;
        for n=1:N1
            a=floor(cc1(n)/delta)+1;
            b=floor(pp1(n)/delta)+1;
            c=floor(dd1(n)/delta)+1;
            if a<=r && b<=r && c<=r
                pdf(a,b,c)=pdf(a,b,c)+1;
            end
        end
        DCP1(:,:,:)=cumsum(cumsum(cumsum(pdf,1),2),3)/N1;

        pdf=zeros(r,r,r);
        cc2(cc2==1)=1-delta/2; pp2(pp2==1)=1-delta/2; dd2(dd2==1)=1-delta/2;
        for n=1:N2
            a=floor(cc2(n)/delta)+1;
            b=floor(pp2(n)/delta)+1;
            c=floor(dd2(n)/delta)+1;
            if a<=r && b<=r && c<=r
                pdf(a,b,c)=pdf(a,b,c)+1;
            end
        end
        DCP2(:,:,:)=cumsum(cumsum(cumsum(pdf,1),2),3)/N2;

        distance=sqrt(sum(sum(sum((DCP1-DCP2).^2)))); 

end

end