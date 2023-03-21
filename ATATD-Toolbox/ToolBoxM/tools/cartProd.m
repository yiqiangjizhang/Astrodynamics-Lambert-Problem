function index = cartProd(n)

% cartProd.m - Cartesian product matrix.
%
% PROTOTYPE:
%   index = cartProd(n)
%
% DESCRIPTION:
%   Computes the cartesian product matrix. Given a number of elements m,
%   and a set of parameters for each element n(m), the cartesian product
%   matrix is the matrix in which each row represents a combination that is
%   possible to make picking, for each element, a parameter from the
%   corresponding set.
%
% INPUT:
%	n[m]                Integer vector with the number of indexes for each 
%                       variable.
%
% OUTPUT:
%   index[prod(n),m]    Cartesian product matrix. One element for each row.
%                       It has a number of columns equal to length(n), and
%                       a number of rows equal to prod(n).
%
% EXAMPLE
%   index = cartprod([2 3 2]) gives:
%       index = [1     1     1
%                1     1     2
%                1     2     1
%                1     2     2
%                1     3     1
%                1     3     2
%                2     1     1
%                2     1     2
%                2     2     1
%                2     2     2
%                2     3     1
%                2     3     2]
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 15/06/2006, MATLAB, cartProd.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 15/06/2006, MATLAB, carprod.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   14/05/2007, REVISION: Camilla Colombo
%   30/12/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

npar=length(n);
nindex=1;
for i=1:npar
    nindex=nindex*n(i);
end
index=ones(nindex,npar);
index(1,npar)=0;
indexcount=0;
for i=1:nindex
    index(i,npar)=indexcount+1;
    j=npar;
    while (j>0)
        if(index(i,j)>n(j))
            index(i:nindex,j)=1;
            index(i:nindex,j-1)=index(i:nindex,j-1)+1;
        end
        j=j-1;
    end
    indexcount=index(i,npar);
end

return