clear;
clc;

% -- Elliptic Curve Domain Parameter, setup. --
EC = DomainParam();
EC.p = 159979;
EC.a = 0;
EC.b = 7;
EC.G = Generate_random_basep(EC.a,EC.b,EC.p);
EC.G = generatorG(EC.a,EC.b,EC.p,EC.G);%[69695,34688]
[EC.N,EC.n] = BSGS(EC.a,EC.b,EC.p,EC.G);
EC.h = cofactor(EC.N,EC.n);
% --------------------------------------------

m = 'Trishika'
sK = 10; % shared secret key integer for encrypt and decrypt.
eM = encode(m,EC.a,EC.b,EC.p,sK)

Alice = Entity();
Bob = Entity();
[Alice.k,Alice.Q] = keyGenerate(EC.a,EC.b,EC.p,EC.n,EC.G);
[Bob.k,Bob.Q] = keyGenerate(EC.a,EC.b,EC.p,EC.n,EC.G);

public_key_valid(EC.a,EC.b,EC.p,EC.n,Alice.Q)
public_key_valid(EC.a,EC.b,EC.p,EC.n,Bob.Q)

[C1,C2] = encryptElgamal(EC.a,EC.b,EC.p,EC.n,Bob.Q,EC.G,eM);
Pm = decryptElgamal(EC.a,EC.b,EC.p,Bob.k,C1,C2);
PmD = decode(Pm,sK)

function eM = encode(m,a,b,p,k) 
% One of Koblitz method's (probabilistic method)
% This function encrypts a plain-text message into a sequence of points on
% the elliptic curve given an auxiallry base parameter known by the two
% communicating parties 'k'.

    m = double(m); % converts each letter of the text into its ASCii values 
    
   for( i=1: length(m) ) % 1st runs through each letter
      for( j=1: m(i)*k+k-1 ) % 2nd try values mk+j (i.e mk+1,mk+2...mk+k-1)
         if(tonelli_shanks((m(i)*k+j)^3 + a*(m(i)*k+j) + b, p)>0) % Here it checks wether TS algo returns a sol to y for (mk+j)^3 + a*(mk+j) + b == (mk+j)^2
             eM(i,1) = m(i)*k+j;
             eM(i,2) = tonelli_shanks((m(i)*k+j)^3 + a*(m(i)*k+j) + b, p);
             break;
         end
      end
   end
end

function mD = decode(mE,k)
% One of Koblitz method's
% This function decrypts an array of points (Only x values) from an
% elliptic curve given an auxiallry base parameter known by the two
% communicating parties 'k'.
    
%    for i=1:length(mE)
%        mD(i) = char(mod((mE(i,1)),128));
%    end
    
    for(i=1:length(mE)) % Runs through each x value in the chipertext.
        [mD(i)] = char(mod(floor((mE(i)-1)/k),128)); % each floor( (xi-1)/k ) is converted into ASCii chars.
    end
end

function [Pm] = decryptElgamal(a,b,p,k,C1,C2)
    [M(1,1),M(1,2)] = DandA(k,a,b,p,C1);
    for i=1:length(C2)
        [Pm(i,1),Pm(i,2)] = pointAddition(a,b,p,[C2(i,1),C2(i,2)],[M(1,1),mod((-M(1,2)),p)]);
    end
end

function [C1,C2] = encryptElgamal(a,b,p,n,Q,G,eM)
    r = randi([1,n-1]);
    disp(r);
    [C1(1),C1(2)] = DandA(r,a,b,p,G);
    [rQ(1),rQ(2)] = DandA(r,a,b,p,Q);
    disp(rQ);
    
    for i=1:length(eM)
        [C2(i,1),C2(i,2)] = pointAddition(a,b,p,rQ,[eM(i,1),eM(i,2)]);
    end
end

function logic = public_key_valid(a,b,p,n,Q)
    if isinf(Q(1)) || isinf(Q(2)) || Q(1) > p || Q(2) > p || isInf(a,b,p,Q) == true  
        logic = false;
        return
    elseif DandA(n,a,b,p,Q) == Inf
        logic = true;
        return
    end
end

function [k,Q] = keyGenerate(a,b,p,n,G)

    k = randi([1,n-1]);
    [Q(1),Q(2)] = DandA(k,a,b,p,G);
end

function h = cofactor(N,n)
    h = N/n;
end

function [G] = generatorG(a,b,p,P)
    
    N = BSGS(a,b,p,P);
    n = randi([1,N]);
    h = N/n;
    
    
    while rem(N,n) > 0 && isPrime(n) ~= true
            n = randi([1,N]);
            h = N/n;
    end
    
    disp(n);
            
    [tokenG(1),tokenG(2)] = DandA(h,a,b,p,P);
        
    while isinf(tokenG(1)) 
        tokenG = Generate_random_basep(a,b,p);
        tokenG = generatorG(a,b,p,tokenG);
    end
        
    G = tokenG;
end

function [N,n] = BSGS(a,b,p,P)

    m = randi([ceil(nthroot(p,4)),ceil(nthroot(p,4))]);
    Pj(1,1) = P(1);
    Pj(1,2) = P(2);
    
    for j=2:m 
        [Pj(j,1),Pj(j,2)] = pointAddition(a,b,p,P,[Pj(j-1,1),Pj(j-1,2)]);
    end

    L = 1;
    [Q(1),Q(2)] = DandA(p+1,a,b,p,P); 
    k = 1;
    M = '';
    
    [k2mP(1,1),k2mP(1,2)] = DandA(k*2*m,a,b,p,P);
    [Qplus(1,1),Qplus(1,2)] = pointAddition(a,b,p,Q,k2mP);
    
    while isInf(a,b,p,k2mP) == 1 % it might occur that we pick m where 2*m*P = O, re-choose m incase. 
        m = randi([floor(nthroot(p,4)),p-1]);
        [k2mP(1,1),k2mP(1,2)] = DandA(k*2*m,a,b,p,P);
        [Qplus(1,1),Qplus(1,2)] = pointAddition(a,b,p,Q,k2mP);
    end
    
    while(k>=1)
        k = k+1;
        [Qplus(k,1),Qplus(k,2)] = pointAddition(a,b,p,[k2mP(1,1),k2mP(1,2)],[Qplus(k-1,1),Qplus(k-1,2)]);
        
        for j=1:length(Pj)
            if Pj(j,1) == Qplus(k,1) && isinf(DandA((p+1+(2*m*k)-j),a,b,p,P))
                M=p+1+(2*m*k)-j;
                break
            elseif mod(-Pj(j,1),p) == Qplus(k,1) && isinf(DandA((p+1+(2*m*k)+j),a,b,p,P))
                M=p+1+(2*m*k)+j;
                break
            end
        end
        
        if M > 0 
            break
        end
    end
    
    k = factor(M);
    primeFactors = k(rem(M,k)==0);
    
    if ~factor(M)
        primeFactors = k;
    end
    
    for i=1:length(primeFactors)
        if(isinf(DandA(floor(M/primeFactors(i)),a,b,p,P)))
            M = floor(M/primeFactors(i));
        end
    end
    
    n = M; 
    L = lcm(L,M);
    
    % if 'p' is below 61 bound [p+1-2√p, p+1+2√p] is too large for the prime.  If left, the algo will find a 'N' which (> the true 'N' (the order)) is divisible by 'M' i.e MP = O.
    if p < 61 
        hasseB = p+1-(floor(sqrt(p))):p+1+(floor(sqrt(p)));
    else
        hasseB = p+1-(floor(2*sqrt(p))):p+1+(floor(2*sqrt(p)));
    end
    
    count = 0;
    
    for j=1:length(hasseB)
        if rem(hasseB(j),L) == 0
            count = count + 1;
        end
    end
        
    if count == 1
       for j=1:length(hasseB)
            if rem(hasseB(j),L) == 0
                N = hasseB(j);
                return
            end
       end
    end
    
    if count > 1
        for j=1:length(hasseB)
            if rem(hasseB(j),L) == 0 
                [J] = hasseB(j);
            end
        end
        N = J;
        return
    end
    
    while count == 0
        N = BSGS(a,b,p,P);
    end
end

function [G_pn] = Generate_random_basep(a,b,p)
% This function outputs a random point on a EC specified y^2 mod p = x^3 + 'a'*x
% + b mod p. This is accomplished by picking random x values from [0 - p-1]
% i.e Z_p until an x value is a quadratic residule (the check for euler
% criteria or if point Inf i.e x == 0). We then work out the corressponding
% y value via the tonelli-shanks algorithm

    x = randi([-(p-1)/2,(p-1)/2]); % pick a random element from Z_p
    z = mod(x^3+(a*x)+b,p); % calculate y^2
    
    while(EulerCriteria(z,p) ~= 1 || tonelli_shanks(z,p) == 0 || x==0) % keep choosing random x's till its a quadratic residue to the EC.
        x = randi([-(p-1)/2,(p-1)/2]);
        z = mod(x^3+(a*x)+b,p);
    end
    
    Gx = mod(x,p);
    Gy = tonelli_shanks(z,p); % Caluculate corressponding y value
    
    % below we want to give 50/50 for returning (x,-y) or (x,y)
    r = randi([0,100]);
    if(mod(r,2)==1) 
        Gy = mod(-Gy,p);
    end
    
    [G_pn] = [Gx,Gy];
end

function [x,y] = DandA(n,a,b,p,P) % scalar point multiplication
% This function returns the coordinates for n*P, via the Double & Add algorithm.
% First 'n' is represented in its binary form. Then check the digits of
% this binary representation of 'n' (we ignore the most significant bit as
% this is the case for n=1 >> 1*P = P and just store in nP immediately), if we come across a '1' then double
% the point 'P' (in [nP]) THEN add point 'P' to the stored point stored in [nP]. If
% we come across a '0' then we ONLY double the point.

    binary = dec2bin(n); % decimal to binary number
    nP = [P(1),P(2)]; % Store point 'P'
    
    if(n==0)
        x = 0;
        y = 0;
        return;
    end
    
    if(n==1) % Case for 1P = P
        x = P(1);
        y = P(2);
        return;
    end
            
    for i = 2 : length(binary) % Starting from msb+1
        if binary(i) == '1' % check if '1' do double
           [nP(1),nP(2)] = pointAddition(a,b,p,nP,nP);
           [nP(1),nP(2)] = pointAddition(a,b,p,P,nP); % If no point O then DO points addition
        else
            [nP(1),nP(2)] = pointAddition(a,b,p,nP,nP);
        end
    end
    
    x = nP(1);
    y = nP(2);
end

function [x,y] = pointAddition(a,b,p,P,Q)

    if isInf(a,b,p,P) == 1 && isInf(a,b,p,Q) == 0
        x = Q(1);
        y = Q(2);
        return
    elseif isInf(a,b,p,Q) == 1 && isInf(a,b,p,P) == 0
        x = P(1);
        y = P(2);
        return
    elseif P(1) == Q(1) && P(2) ~= Q(2) || P(1) == Q(1) && P(2) == 0
        x = Inf;
        y = Inf;
        return
    elseif isinf(P(1)) && isinf(Q(1))
        x = Inf;
        y = Inf;
    end
    
    if P(1) ~= Q(1)
       lamda = mod( mod(Q(2) - P(2),p) * modInverse(mod(Q(1) - P(1),p) ,p) ,p); % slope
       x = mod(mod(lamda^2,p) + mod(- P(1) - Q(1),p),p);
       y = mod(-P(2) + mod(lamda * mod(P(1) - x,p),p),p);
       return
    elseif P(1) == Q(1) && P(2) == Q(2) && P(2) ~= 0 && ~isinf(P(1)) && ~isinf(Q(1))
        lamda = mod( (3 * (P(1)^2) + a) * modInverse(2*P(2) ,p) ,p);
        x = mod( mod(lamda^2,p) - (2 * P(1)) ,p);
        y = mod( -P(2) + (lamda * (P(1) - x)) ,p);
        return
    end
end

function inv = modInverse(a,p) 
% This method returns the multiplicative inverse mod p via a recursive method of the extended 
% Euclidian algorithm. Algorithm Assumption: a and p are coprimes, i.e.,
% gcd(a, p) = 1;
    m = p;
    y = 0;
    x = 1;
    
    if(p==1)
        inv = 0;
        return;
    end
    
    while (a > 1)
        q = floor(a / p); % q is quotient.
        t = p;
        
        % p is remainder now, process same as Euclid's algo.
        p = mod(a,p);
        a = t;
        t = y;
        
        % Update x and y.
        y = x - q * y;
        x = t;
    end
    
    % make x positive.
    if(x<0)
        x = x + m;
    end
    inv = x;
end

function q = powerMod(b,e,m)
% This function is the fast expoentiation which outputs = b^e mod m using the binary modular expoentiation
% algorithm. If there is a '1' in the bin number we sqaure and multiply, if
% '0' we sqaure.

    binary = dec2bin(e);
    result = mod(b,m); % Since the 1st binary digit is 1, simply list the number. 
    
    for (i = 2 : length(binary))
        if(binary(i) == '1')
            result = mod((result^2), m);
            result = mod(result * b, m);
        else
            result = mod(result^2, m);

        end
    end
    q = result;
end

function rs = EulerCriteria(a,p)
% This function outputs logical values for the outcome of Euler's
% Crtierion: a^((p-1)/2) === 1 mod p. i.e. 'a' is a quadratic residual of
% 'p' (odd prime specifically). If Zero(mod p) is passed then return true 
% since 0 is always a quadratic residual. We don't need to check if
% gcd(a,p) == 1;

    base = a;
    power = ((p-1)/2);
    
    if (powerMod(base,power,p) == 1 || mod(a,p) == 0)
        rs = true;
    else
        rs = false;
    end
end

function y = tonelli_shanks(n,p)
% This function finds/returns a mod-p square root of 'n', given an odd 
% prime 'p' and integer 'n' i.e. y = sqrt(n) mod p. Assume n is coprime to p gcd
% = 1.

    % Case when p|n, so n = 0 (mod p). 
    if(mod(n,p) == 0)
        y = 0;
        return;
    end
    if(EulerCriteria(n,p) == 0)
        y = -1;
        return;
    end
    % case if p = 3 mod 4.    
    if (mod(p,4) == 3) 
        y = powerMod(n, ((p+1)/4), p);
        return;
    end
    % case if p = 1 mod 4.
    Q = p-1;
    S = 0;
    % factor out powers of 2 to find Q & S s.t. p-1 = Q2^S where 'Q' is odd
    while (mod(Q,2) == 0)
        S = S+1;
        Q = Q/2;
    end
    
    % find non-residual of p by brute.
    z = randi([-(p-1)/2,(p-1)/2]);
    while(EulerCriteria(z,p) == 1)
        z = randi([-(p-1)/2,(p-1)/2]);
    end
    
    M = S;
    c = powerMod(z, Q, p);
    t = powerMod(n, Q, p);
    R = powerMod(n, floor((Q+1)/2), p);
    
    while(t ~= 1)
        i = 0;
        temp = t;
        while(temp ~= 1)
            i = i+1;
            temp = mod(temp*temp,p);
        end
        pow2 = 2^(M-i-1);
        b = powerMod(c, pow2, p);
        M = i;
        c = mod(b*b,p);
        t = mod(t * b*b,p);
        R = mod(R * b, p);
    end
    y = R;
    return;
end

function probp = isPrime(n)
% Fermat's Little Theorem.
% A probabilistic test to determine whether a number is a probable prime.
% If 'p' is prime and 'a' is not divisible by 'p', then 'a^(p-1) mod p = 1'.
% We pick random integers 'a' not divisible by 'p' and see whether the equality holds.

    if (n <= 1 || n == 4) 
        probp = false;
        return;
    end
    if (n <= 3) 
        probp = true;
        return;
    end
    
    k = floor(sqrt(n));
    % Try k times
    while (k>0)
        % Pick a random number in [2..n-2]
        a = randi([2,n-2]); 
        % Fermat's little theorem
        if(powerMod(a,n-1,n)~=1)
            probp = false;
            return;
        end
        k = k-1;
    end
    probp = true;
    return;
end


function logic = isInf(a,b,p,P)
    if mod(P(1)^3 + a*P(1) + b,p) == mod(P(2)^2,p)  % check is final point calculated is the point O return Inf
        logic = false;
        return
    else
        logic = true;
        return
    end
end