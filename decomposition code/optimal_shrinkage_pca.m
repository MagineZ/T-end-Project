function [eigenvals sigma] = optimal_shrinkage(eigenvals,gamma,loss,sigma)

% function singvals = optimal_shrinkage(singvals,gamma,sigma_known)
%
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% Perform optimal shrinkage (w.r.t one of a few possible losses) on data
% eigenvalues, when the noise is assumed white, and the noise level is known 
% or unknown. Use analytical formulas for optimal shrinkage.
%
% IN:
%   eigenvals: a vector of data eigenvalues, obtained by running eig
%              on the p-by-p sample covariance matrix corresponding to a dataset
%              of n samples in dimension p 
%
%   gamma:     aspect ratio p/n of the dataset. For convenience we assume p<=n.
%             
%   loss:     the loss function for which the shrinkage should be optimal
%             (see paper). 
%
%             presently implemented: 
%             F_1: Frobenius norm on A-B
%             F_2: Frobenius norm on precision (A^-1 - B^-1)
%             F_3: Frobenius norm on A^-1 * B - I
%             F_4: Frobenius norm on B^-1 * A - I
%             F_6 Frobenius norm on A^1/2 * B * A^1/2 - I
%             N_1: Nuclear norm on A-B
%             N_2: Nuclear norm on precision (A^-1 - B^-1)
%             N_3: Nuclear norm on A^-1 * B - I
%             N_4: Nuclear norm on B^-1 * A - I
%             N_6 Nuclear norm on A^1/2 * B * A^1/2 - I
%             O_1: Operator norm on A-B
%             O_2: Operator norm on precision (A^-1 - B^-1)
%             O_6 Operator norm on A^1/2 * B * A^1/2 - I
%             Stein: Stein's loss
%             Ent: Entropy loss
%             Div: Divergence loss
%             Fre: Frechet loss
%             Aff: Affine loss
%
%   sigma:    (optional) noise standard deviation (of each entry of the noise matrix) 
%             if this argument is not provided, the noise level is estimated 
%             from the data.
%
% OUT:
%    eigenvals: the vector of eigenvalues after performing optimal shrinkage
%    sigma: an estimate of the noise level
%
% Usage:
%   Given an n-by-p data matrix Y assumed to follow approximately the spiked
%   model (that is, each row is an i.i.d sample from a p-variate Gaussian
%   distribution, whose population covariance is a multiple of the
%   identity except for a small number of top population "spikes", we form an
%   estimate of the poplation covariance as follows -
%
%   [n,p] = size(Y);
%   S = Y'*Y;
%   [V D] = eig(S);
%   d = diag(D);
%   d = optimal_shrinkage(d , p/n ,'F_1');  
%   Shat = V * diag(d) * V';
%
%   where you can replace 'F_1' with one of the other losses. 
%   if the noise level sigma is known, in the fifth line use instead e.g
%       d = optimal_shrinkage( d , p/n , 'F_1' , sigma);
%    
% -----------------------------------------------------------------------------
% Authors: David Donoho and Matan Gavish <lastname>@stanford.edu, 2015
% 
% This program is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------

assert(prod(size(gamma))==1)
assert(gamma<=1)
assert(gamma>0)
assert(prod(size(eigenvals))==length(eigenvals))

% estimate sigma if needed
if nargin<4
    warning('off','MATLAB:quadl:MinStepSize')
    MPmedian = MedianMarcenkoPastur(gamma);
    sigma = sqrt(median(eigenvals) / MPmedian);
    fprintf('estimated sigma=%0.2f \n',sigma);
end

eigenvals = optshrink_impl(eigenvals,gamma,loss,sigma);

end


function eigenvals = optshrink_impl(eigenvals,gamma,loss,sigma)

    sigma2 = sigma^2;
    lam_plus = (1+sqrt(gamma))^2;
    
    ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) + sqrt((lam+1-gamma).^2-4*lam))/(2.0));

    c = @(lam)((lam>=lam_plus) .* sqrt( (1-gamma./((ell(lam)-1).^2) ) ./ (1+gamma./(ell(lam)-1) )   ));
    s = @(lam)(sqrt(1-c(lam).^2 ) );

    impl_F_1 = @(ell,c,s)(max(1+(c.^2).*(ell-1),0));
    impl_F_2 = @(ell,c,s)(max(ell./((c.^2)+ell.*(s.^2)),0));
    impl_F_3 = @(ell,c,s)(max(1 + (ell-1).*((c.^2)./((ell.^2).*(s.^2) + c.^2 )),1));
    impl_F_4 = @(ell,c,s)((s.^2 + (ell.^2).*(c.^2))./((s.^2) + (ell.*(c.^2))));
    impl_F_6 = @(ell,c,s)(1 + ((ell-1).*(c.^2)) ./ (((c.^2)+ell.*(s.^2)).^2));
    impl_O_1 = @(ell,c,s)(ell); % this is just debiasing to pop eigen
    impl_O_2 = @(ell,c,s)(ell); % this is just debiasing to pop eigen
    impl_O_6 = @(ell,c,s)(1+((ell-1)./(c.^2 + ell.*(s.^2)))); 
    impl_N_1 = @(ell,c,s)(max(1+(ell-1).*(1-2*(s.^2)),1));
    impl_N_2 = @(ell,c,s)(max(ell./((2*ell-1).*(s.^2) + c.^2),1));
    impl_N_3 = @(ell,c,s)(max(ell./(c.^2+(ell.^2).*(s.^2)),1));
    impl_N_4 = @(ell,c,s)(max(ell .* (c.^2) + (s.^2) ./ ell,1));
    impl_N_6 = @(ell,c,s)(max((ell - ((ell-1).^2) .* ...
                           (c.^2) .* (s.^2))./(((c.^2) + ell.*(s.^2)).^2),1));
    impl_Stein =  @(ell,c,s)(ell./(c.^2 + ell.*(s.^2))); 
    impl_Ent =  @(ell,c,s)(ell.*(c.^2) + s.^2); 
    impl_Div =  @(ell,c,s)(sqrt(((ell.^2).*(c.^2)+ell.*(s.^2))./(c.^2+(s.^2).*ell))); 
    impl_Fre =  @(ell,c,s)((sqrt(ell).*(c.^2)+s.^2).^2); 
    impl_Aff =  @(ell,c,s)( ((1+c.^2).*ell + (s.^2)) ./ (1+(c.^2)+ell.*(s.^2))); 

    assert(sigma>0)
    assert(prod(size(sigma))==1)
    eigenvals = eigenvals / sigma2;
    I = (eigenvals > lam_plus);
    eigenvals(~I) = 1;
    str1 = [loss '_func = @(lam)(max(1,impl_' loss '(ell(lam),c(lam),s(lam))));'];
    str2 = ['eigenvals(I)= ' loss '_func(eigenvals(I));'];
    eval(str1);
    eval(str2);
    eigenvals = sigma2 * eigenvals;
end


function I = MarcenkoPasturIntegral(x,gamma)
    if gamma <= 0 | gamma > 1,
        error('gamma beyond')
    end
    lobnd = (1 - sqrt(gamma))^2;
    hibnd = (1 + sqrt(gamma))^2;
    if (x < lobnd) | (x > hibnd),
        error('x beyond')
    end
    dens = @(t) sqrt((hibnd-t).*(t-lobnd))./(2*pi*gamma.*t);
    I = quadl(dens,lobnd,x);
    fprintf('x=%.3f,gamma=%.3f,I=%.3f\n',x,gamma,I);
end


function med = MedianMarcenkoPastur(gamma)
    MarPas = @(x) 1-incMarPas(x,gamma,0);
    lobnd = (1 - sqrt(gamma))^2;
    hibnd = (1 + sqrt(gamma))^2;
    change = 1;
    while change & (hibnd - lobnd > .001),
      change = 0;
      x = linspace(lobnd,hibnd,5);
      for i=1:length(x),
          y(i) = MarPas(x(i));
      end
      if any(y < 0.5),
         lobnd = max(x(y < 0.5));
         change = 1;
      end
      if any(y > 0.5),
         hibnd = min(x(y > 0.5));
         change = 1;
      end
    end
    med = (hibnd+lobnd)./2;
end

function I = incMarPas(x0,gamma,alpha)
    if gamma > 1,
        error('gammaBeyond');
    end
    topSpec = (1 + sqrt(gamma))^2;
    botSpec = (1 - sqrt(gamma))^2;
    MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
                         sqrt((topSpec-x).*(x-botSpec))./(gamma.* x)./(2 .* pi), ...
                         0);
    if alpha ~= 0,
       fun = @(x) (x.^alpha .* MarPas(x));
    else
       fun = @(x) MarPas(x);
    end
    I = quadl(fun,x0,topSpec);
    
    function y=IfElse(Q,point,counterPoint)
        y = point;
        if any(~Q),
            if length(counterPoint) == 1,
                counterPoint = ones(size(Q)).*counterPoint;
            end
            y(~Q) = counterPoint(~Q);
        end
        
    end
end


