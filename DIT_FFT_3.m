close all ; clear all ; clc ;

% Taking the input from the user ,[1 2 3....]
x=input('Enter the sequence in []: ');    

% The length of the input x must be an integer power of 2. 
x_rounded=nextpow2(length(x));  

x_new=[x zeros(1,(2^x_rounded)-length(x))];   % Zero padding        
N=length(x_new);                              % Length of zero padded sequence                           
S=log2(N);                                    % Number of stages                          
% The input to the butterfly structure is in bitreverse form of the input 
x_new=bitrevorder(x_new);

% Let's take an example of 8-point 
% There will be 3 stages .i.e 1st stage-2point dft : 2nd-4pt : 3rd-8pt 
for stage=1:S
   p=1;                                    % p -> start point for 1st part
   q=1+2^(stage-1);                        % q -> start point for 2nd part
   n=0;                                    % n -> Twiddle factor for each iteration
   
   % p & q are used as indexes
   % This while loop is used for each iteration in a stage
   % ex: in stage 1 loop is repeated 4 times , in stage 2 loop is repeated 2 times and so on 
   while( n<=2^(stage-1)-1 && q<=N)
       e=exp((-1i)*2*pi*n/(2^stage));      % Twiddle factor calculation
       G=x_new(p)+e*x_new(q);              % G and H are calculated , which are the output of the current stage   
       H=x_new(p)-e*x_new(q);
       x_new(p)=G;                         % G and H are used as input for the next stage 
       x_new(q)=H;
       p=p+1;                              % Incrementing p,q & n for the next iteration of same stage
       q=q+1;     
       n=n+1;
       if(rem(q,2^stage)==1)               % Once all the computations of the iteration are completed
          p=p+2^(stage-1);                 % Incrementing p,q and resetting n for the next butterfly structure 
          q=q+2^(stage-1);
          n=0;
       end
   end
   % x_new will be printed as many times as number of stages 
   % and will be in order from x(n) till we reach X(k)
   x_new
end

X=fft(x,N)  % calculating fft using built-in function for verification at last stage
