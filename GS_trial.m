function gamma=GS_trial(startP,endP)
a=startP;
b=endP;
GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
d=GR*(b-a)+a;
c=b-GR*(b-a);
% ff=cell(1,5);%prelocation for speed
ff=@(x)x.^2+2*x+3;

while abs(c-d)>0.00001
    %     for i=1:5
    %         ffc=ffc+ff{i}(c);
    %         ffd=ffd+ff{i}(d);
    %     end
    ffc=ff(c);
    ffd=ff(d);
    if ffc<=ffd
        b=d;
        d=c;
        c=b-GR*(b-a);
    else
        a=c;
        c=d;
        d=a+GR*(b-a);
    end
end
gamma=(a+b)/2;
end