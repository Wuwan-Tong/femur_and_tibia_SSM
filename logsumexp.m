function s = logsumexp(x, dim)
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).
% if nargin == 1, 
%     % Determine which dimension sum will use
%     dim = find(size(x)~=1,1);
%     if isempty(dim), dim = 1; end
% end

% subtract the largest in each row
y = max(x,[],1);
x = bsxfun(@minus,x,y);
s = log(sum(exp(x),dim));
s = bsxfun(@plus,y,s);
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end
