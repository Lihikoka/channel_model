function E=erfcomp(z);  
% ERFCOMP	(komplexes) Fehlerintegral  
%  
%		E=ERFCOMP(Z)  
%  
%		Die (komplexe) Matrix Z ist das Eingabeargument  
%		f?r diese Prozedur. Zur?ckgeliefert wird die  
%		komplexe Error-Funktion.  
%  
%		Entnommen aus Abramowitz,  
%		Abschnitt "Error Functions and Fresnel Integrals",  
%		Seite 299  
  
[zx,sx]=size(z);				% Dimension der Eingabematrix  
x=real(z); y=imag(z);  
x=x(:); y=y(:);  
  
i=sqrt(-1);  
						% E=erf(x);  
KOa=+0.3275911;					% Hier erfolgt die Berechnung der  
KOA=+0.254829592;				% Error-Funktion durch eine rationale  
KOB=-0.284496736;				% Approximation f?r reelle Argumente.  
KOC=+1.421413741;  
KOD=-1.453152027;  
KOE=+1.061405429;  
w=1 ./(1 + KOa .*abs(x));  
P=KOE .*w + KOD;  
P=P .*w + KOC;  
P=P .*w + KOB;  
P=P .*w + KOA;  
P=P .*w;  
Q=exp(-(x .*x));  
E=(1 - Q .*P) .*sign(x);  
  
S=zeros(size(E));  
abserr=100;  
n=1;  
while abserr>1E-12,  
	fn=2 .*x - 2 .*x .*cosh(n .*y) .*cos(2 .*x .*y) + n .*sinh(n .*y) .*sin(2 .*x .*y);  
	gn=2 .*x .*cosh(n .*y) .*sin(2 .*x .*y) + n .*sinh(n .*y) .*cos(2 .*x .*y);  
	fgn=(fn + i .*gn) .*exp(-0.25 .*(n .^2)) ./( (n .^2) + (4 .*(x .^2)));  
	Sakt=2 .*fgn .*exp(-(x .^2)) ./pi;  
	abserr=max(abs(Sakt));  
	n=n+1;  
	S=S + Sakt;  
end  
  
ki=find(x==0 & y==0);  
if ~isempty(ki),  
	E(ki)=zeros(size(ki));  
end  
ki=find(x==0 & y~=0);  
if ~isempty(ki),  
	E(ki)=i .*y(ki) ./pi;  
end  
ki=find(x~=0 & y~=0);  
if ~isempty(ki),  
	xk=x(ki); yk=y(ki);  
	E(ki)=erf(xk) + (exp(-(xk .^2)) .*( (1 - cos(2 .*xk .*yk)) + i .*sin(2 .*xk .*yk)) ./(2 .*pi .*xk));  
end  
E=E+S;  
E=reshape(E,zx,sx);
