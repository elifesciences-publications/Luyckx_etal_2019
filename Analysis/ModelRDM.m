function mods = ModelRDM
%(c) Bernhard Spitzer, 2016

%% Numerical distance

pc=zeros(6,6);
for i=1:6
    for j=1:6
        pc(i,j)=abs(i-j)./5;
    end
end

mods.numd=pc;

end