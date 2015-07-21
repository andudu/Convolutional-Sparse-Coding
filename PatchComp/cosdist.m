function [ a] = cosdist( y,z )
%COSDIST Summary of this function goes here
%   Detailed explanation goes here

if size(y,1) == 1
    y = y';
end

if size(z,1) == 1
    z = z';
end


a = (y'*z)/(norm(z)*norm(y));
a = abs(a);
end

