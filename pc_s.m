function [s_in]=pc_s(pc_in,s_im)
s=s_im:0.001:0.9;
pc=-0.00011*exp(-44.02*(s-0.496))+278.3*exp(8.103*(s-0.496))-191.8;

for i=1:length(pc_in)
    if pc_in<min(pc)
        l=max(2,find(pc>=pc_in(i),1));
        s_in(i)=max(s_im,(s(l)+s(l-1))/2);
    else
    s_in(i)=s_im;
    end
end
end
