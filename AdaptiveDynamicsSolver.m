% Given the mutation variance covariance matrix and the selection gradient,
% which we assume can be calculated analytically, 
% numerically integrate the adaptive dynamics


M = eye(d);
s[i] = sum(b[i,j]*x[j],d)+sum(sum(a[i,j,k]*x[j]*x[k],j),k);