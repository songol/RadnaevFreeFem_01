/*Dir-Neum*/
/*Radnaev T.Ts.*/
real l = 1.;
int nb = 100;
real[int] yout = [sin(pi/3), sin(2*pi/3), sin(pi), sin(4*pi/3),sin(5*pi/3), sin(2*pi)] * l;
real[int] xout = [cos(pi/3), cos(2*pi/3), cos(pi), cos(4*pi/3),cos(5*pi/3), cos(2*pi)] * l;

real[int] yin = [sin(0.0), sin(2*pi/3), sin(4*pi/3)] * (l / 2.0);
real[int] xin = [cos(0.0), cos(2*pi/3), cos(4*pi/3)] * (l / 2.0);

int[int] nout = [1, 1, 1, 1, 1, 1] * nb;
int[int] nin = [1, 1, 1] * nb;

int[int] labelout = [0, 1, 2, 3, 4, 5];
int[int] labelin = [6, 7, 8];

border GammaOut(t = 0.0, 1.0; i) {
	int ii = (i + 1)%6;
	real t1 = 1 - t;
	x = xout[i]*t1 + xout[ii]*t;
	y = yout[i]*t1 + yout[ii]*t;
	label = labelout[i];
}

border GammaIn(t = 1.0, 0.0; i){
	int ii = (i + 1)%3;
	real t1 = -1 + t;
	x = xin[i]*t1 - xin[ii]*t;
	y = yin[i]*t1 - yin[ii]*t;
	label = labelin[i];
}

//plot(GammaOut(nout) + GammaIn(nin));

mesh Th1 = buildmesh(GammaOut(nout) + GammaIn(nin));
//plot(Th1, ps = "Mesh");

fespace Vh(Th1, P1);
Vh u, v;
Vh uexact = sin(3*x - 2*y); 

/*g = (du/dn) on top and bottom parts of outer border */

//func gtop = -2*cos(3*x - 2*y);
//func gbottom = 2*cos(3*x - 2*y);
func g = 3*cos(3*x - 2*y) * N.x - 2*cos(3*x - 2*y) * N.y;


func f = 13*sin(3*x - 2*y);
solve solution(u, v) = int2d(Th1)(dx(u)*dx(v) + dy(u)*dy(v))
					- int2d(Th1)(f*v)
					//- int1d(Th1, labelout[3])(gbottom*v)
					//- int1d(Th1, labelout[0])(gtop*v)
					-int1d(Th1, labelout[3])(g*v)
					-int1d(Th1, labelout[0])(g*v)
					+ on(labelout[1], u = uexact)
					+ on(labelout[2], u = uexact)
					+ on(labelout[4], u = uexact)
					+ on(labelout[5], u = uexact)
					+ on(labelin[2], u = uexact)
					+ on(labelin[1], u = uexact)
					+ on(labelin[0], u = uexact)
					;
real[int] pogr(6);
//Vh uerror = abs(u - uexact);
//real uerrorL2 = sqrt(int2d(Th1)(uerror^2));
for (int i = 0; i < 6; i++)
{
	int nb = 5 * 2^i;
	real l = 1.;
	int[int] nout = [1, 1, 1, 1, 1, 1] * nb;
	int[int] nin = [1, 1, 1] * nb;
	mesh Th1 = buildmesh(GammaOut(nout) + GammaIn(nin));
	fespace Vh1(Th1, P1);
	Vh1 u, v;
	Vh1 uexact = sin(3*x - 2*y); 

	solve solution(u, v) = int2d(Th1)(dx(u)*dx(v) + dy(u)*dy(v))
					- int2d(Th1)(f*v)
					//- int1d(Th1, labelout[3])(gbottom*v)
					//- int1d(Th1, labelout[0])(gtop*v)
					-int1d(Th1, labelout[3])(g*v)
					-int1d(Th1, labelout[0])(g*v)
					+ on(labelout[1], u = uexact)
					+ on(labelout[2], u = uexact)
					+ on(labelout[4], u = uexact)
					+ on(labelout[5], u = uexact)
					+ on(labelin[2], u = uexact)
					+ on(labelin[1], u = uexact)
					+ on(labelin[0], u = uexact)
					;
	
	Vh uerror = abs(u - uexact);
	real uerrorL2 = sqrt(int2d(Th1)(uerror^2));
	pogr[i] = uerrorL2;
}

for (int i = 0; i < 6; i++){
	cout << "error = " << pogr[i] << "; N = " << 5 * 2^i << endl;
}
//plot(u, ps = "Numeric_solution", fill = 1, wait = 1, value = 0.5, nbiso = nb);
//plot(uexact, ps = "Exact_solution", fill = 1, value = 0.5, nbiso = nb);
//plot(uerror, ps = "Abs_error", fill = 1, value = 0.01, nbiso = nb);

  