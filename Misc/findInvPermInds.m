function invPermInds = findInvPermInds(permInds)
% function to find the indices that would reverse a permutation

permTable = sortrows([permInds' (1:length(permInds))']);
invPermInds = permTable(:,2);