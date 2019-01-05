 function onIT = onLINE(v,A,B)
%% onLINE finds out if v=(x,y) is on the line segment AB
% This function requires A and B to have 2 elements. Further to this, it
% allows v to be (2,N) or (N,2) in size -- where N is a positive integer.
% This function returns a vector of length N which stores true values
% wherever v is on/in AB.


tol=1.0e-06; % the relative tolerance to use in deciding if v is on or off 
             % the line segment AB 

if (numel(A) ~= 2) || (numel(B) ~= 2)  %make sure only two elements in A and B
    disp('Either A does not have 2 elements, or B does not')
    onIT=false;
    return
end

if size(A,1) == 1   % make sure A is a column vector. Presently it has only one row.
    A=A';
end

if size(B,1) == 1   % make sure B is a column vector
    B=B';
end


% dot treats columns as vectors to be dotted:
if size(v,1) ~= 2    % watch out for v not having two rows
    if size(v,2) == 2  % it has two columns, so tranpose it
        v=v';
    else
          disp('v does not have 2 rows, or 2 cols, so bailing from onLINE')
          onIT=false; display(size(v));
          return
    end
end 

% project v onto the vector AB:
    
vAB=B-A;  % the vector from point A to point B
magAB=norm(vAB);  % the magnitude of the vector AB
tol = tol*magAB; % the absolute tolerance for deciding


N=size(v,2); % number of columns of v. It must have two rows.
arrA=repmat(A,1,N);
arrVmA= v-arrA; arrvAB=repmat(vAB,1,N);

% dot will take the dot product of corresponding column vectors and produce
% a row vector:
dotpVec = dot(arrVmA,arrvAB); % a row vector of the dot products
dotparr= repmat(dotpVec,2,1); % produce a second row of the same values dot prod

Projection = (arrvAB .* dotparr)/(magAB*magAB);  % the vector v projected onto AB (each column)
Orthog = arrVmA - Projection; % the vector from v, orthogonal to the line through AB


Checknorm = hypot(Orthog(1,:), Orthog(2,:)) > tol; % true if the orthogonal distance is greater than the tolerance so not on AB
                                % a row vector of 0, 1 or false, true

onIT = true(1,N); % set up true for every column vector v initially
onIT(Checknorm) = false;      % set to false, all cols with norms > tol

t = (v(1,:) - arrA(1,:))./arrvAB(1,:); % either component can be used here, BTW. t is a row vec, N long


 indexOFF = or(t < -tol/magAB, t>1+tol/magAB);  % off the line if t is outside range [0,1]
 onIT(indexOFF)=false; % set any places in line but with t outside (0,1) to false; not on AB


end