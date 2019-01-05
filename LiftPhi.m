function theta=LiftPhi(phi, H, phi0, Hstar, A,B,C,D)
%% LiftPhi converts phi into theta and returns the associated theta value.

% Note: that each phi value has two components, phi1 and phi2 such that
% phi=[phi1; phi2]; H is the height of the bauxite; phi0 is the vector value
% of phi when there is no bauxite present; Hstar is 1x3 and contains the
% values of bauxite height that allow to compute which Riemann surface the 
% phase is on: Hstar=[H_1^*, H_2^*, H_{2pi}] (with dimensions 2x1); A,B,C,
% and D are the vertices of the parallelogram. 
% in general it is a vector with one row but as many columns as phi. Theta is 
% the lift of phi.

% If the order of phi1 and phi2 is reversed, this function will return the
% same value for theta. H and phi can have more than one column if one wishes
% to calculate multiple theta values simultaneously, in this case, theta
% will have two rows, and the same number of columns as H and phi. 

     H1s=Hstar(1); H2s = Hstar(2); H2pi=Hstar(3);
     H1ss = H1s+H2pi/2; H2ss = H2s+H2pi/2;
     
     k1=ceil((H-H1s)/H2pi); n1= ceil((H-H1ss)/H2pi);  
     k2=ceil((H-H2s)/H2pi); n2= ceil((H-H2ss)/H2pi);  % this also is a vector, 1xlength(H)
     PhiInAB = onLINE(phi,A,B);  PhiInBC = onLINE(phi,B,C); 
     PhiInCD = onLINE(phi,C,D);  PhiInDA = onLINE(phi, D, A);
     theta = zeros(size(phi));
      
     P1ind=PhiInAB|PhiInBC; P2ind=PhiInCD | PhiInDA;
     
     
     theta(1,P1ind) = phi0(1) - phi(1,P1ind) + 2*k1(P1ind)'*pi;
     theta(1,P2ind) = phi0(1) + phi(1,P2ind) + 2*n1(P2ind)'*pi;

     P3ind= PhiInCD | PhiInBC; P4ind=PhiInAB | PhiInDA;
     theta(2,P3ind) = phi0(2) - phi(2,P3ind) + 2*k2(P3ind)'*pi;
     theta(2,P4ind) = phi0(2) + phi(2,P4ind) + 2*n2(P4ind)'*pi;
     
end