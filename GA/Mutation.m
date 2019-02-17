function mut_Chil = Mutation(parents,~,~,~,~,~,this_Pop,Constrs)
% Custom Mutation function for G.A. search
PTB = rand(size(this_Pop(parents,:)))>0.8;
while min(max(PTB,[],2))==0
    J = find(max(PTB,[],2)==0);
    PTB(J,:) = rand(size(PTB(J,:)))>0.8;
end

mut_Chil = this_Pop(parents,:)+...
    PTB.*randn(size(this_Pop(parents,:)))/10^(5*rand)+...
    10^(-3*rand)*randn*(rand(size(this_Pop(parents,:)))>0.95);

% This mutation function takes the original parents, then chooses 50% of
% the number to mutate, then mutates these by a normally distributed random
% variable (multiplicatively). Finally, we add another small normally distributed
% random variable (again to 50%) in order to push the values away from
% zero.

FLP = ones(size(mut_Chil))-2*(rand(size(mut_Chil))>0.99);
mut_Chil=mut_Chil.*FLP;
for i=1:size(mut_Chil,1)
    mut_Chil(i,:) = max([mut_Chil(i,:);Constrs.LB]);
    mut_Chil(i,:) = min([mut_Chil(i,:);Constrs.UB]);
end

end

