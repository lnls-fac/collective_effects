function frac_cont(num, n)

c = num;
a = zeros(1,n);
b = zeros(1,n);
for i=1:n
    a(i) = mod(c,1);
    b(i) = c-a(i);
    c= a(i)^-1;
end

num3 = b(n);
for i=n-1:-1:1
    num3 = 1/num3 + b(i);
end


fprintf('%4d',b);
fprintf('\n%7.5g\n',num-num3);
