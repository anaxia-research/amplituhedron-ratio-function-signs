"""Symbolic checks for the y-parametrization of Region I."""
import sympy as sp
a,b,c=sp.symbols('a b c',positive=True)
u=a*(1-b)*(1-c)/((1-a*b)*(1-a*c)); v=b*(1-c)*(1-a)/((1-b*c)*(1-b*a)); w=c*(1-a)*(1-b)/((1-c*a)*(1-c*b))
D=(1-a*b)*(1-b*c)*(1-c*a)
s=u+v+w-1; P=u*v*w; r=(1-a)*(1-b)*(1-c)*(1-a*b*c)/D
chk=lambda e: sp.simplify(sp.together(e))==0
print('s formula', chk(s+(1-a)*(1-b)*(1-c)*(1+a*b*c)/D))
print('Delta = r^2', chk(s**2-4*P-r**2))
zp=(s+r)/2; zm=(s-r)/2
print('z+ formula',chk(zp+a*b*c*(1-a)*(1-b)*(1-c)/D),'z- formula',chk(zm+(1-a)*(1-b)*(1-c)/D))
print('inverse a',chk((u-zp)/(u-zm)-a),'b',chk((v-zp)/(v-zm)-b),'c',chk((w-zp)/(w-zm)-c))
print('-x_u^+',chk(-u*zp/P-a*(1-b*c)/(1-a)),' -x_u^-',chk(-u*zm/P-(1-b*c)/(b*c*(1-a))))
print('x+/x-',chk(zp/zm-a*b*c),' 1/u-1',chk(1/u-1-(1-a)*(1-a*b*c)/(a*(1-b)*(1-c))))
print('(1-u)(1-v)(1-w) identity', chk((u*v+v*w+w*u)-s-P-(1-u)*(1-v)*(1-w)))
# converse: Phi(psi(u,v,w)) = (u,v,w) using r^2=Delta
U,V,W,R=sp.symbols('U V W R',positive=True)
S=U+V+W-1; PP=U*V*W
ZP=(S+R)/2; ZM=(S-R)/2
A=(U-ZP)/(U-ZM); B=(V-ZP)/(V-ZM); Cc=(W-ZP)/(W-ZM)
uu=A*(1-B)*(1-Cc)/((1-A*B)*(1-A*Cc))
num=sp.numer(sp.together(uu-U))
rem=sp.rem(sp.Poly(sp.expand(num),R).as_expr(), R**2-(S**2-4*PP), R)
print('Phi(psi)=id (u-component):', sp.expand(rem)==0)
