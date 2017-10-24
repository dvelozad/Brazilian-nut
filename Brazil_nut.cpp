#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
#define M_PI 3.14159265358979323846
using namespace std;
const double g = 9.80;
const double K = 1e4, Gamma = 50, Kcundall = 10, MU = .4, k_Shake = 1e2, Amplitud = -2*6/175;

//Boundary dimensions - 0 cubical box
//				      - 1 cilindrical box
const double Lx = 30, Ly = 30, Lz = 175;
const double R_cil = 24;
const int boundary = 1;


const int Nx = 1, Ny = 1, Nz = 28;
const int Nr = 6, Ntheta = 30, N = Nr*Ntheta*Nz;

const double Zeta = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

class Cuerpo;
class Colisionador;

//----------------------Clase Cuerpo---------------
class Cuerpo {
private:
	vector3D r, V, F, omega, tau; double I, m, R, theta;
public:
	void Inicie(double x0, double y0, double z0,
		double Vx0, double Vy0, double Vz0,
		double theta0,
		double omegax0, double omegay0, double omegaz0,
		double m0, double R0);
	void BorreFuerza_Torque(void);
	void AgregueFuerza(vector3D F0);
	void AgregueTorque(vector3D tau0);
	void Mueva_r(double dt, double Constante);
	void Mueva_V(double dt, double Constante);
	void Dibujese(void);
	double Getx(void) { return r.x(); }; //Función Inline
	double Gety(void) { return r.y(); }; //Función Inline
	double Getz(void) { return r.z(); }; //Función Inline
	double Get_r(void) { return norma(r); }; //Función Inline
	double Get_R(void) { return R; }; //Función Inline
	friend class Colisionador;
};
void Cuerpo::Inicie(double x0, double y0, double z0,
	double Vx0, double Vy0, double Vz0,
	double theta0,
	double omegax0, double omegay0, double omegaz0,
	double m0, double R0) {
	r.cargue(x0, y0, z0); V.cargue(Vx0, Vy0, Vz0); m = m0; R = R0;  omega.cargue(omegax0, omegay0, omegaz0);
	theta = theta0;  m = m0;  R = R0;  I = (2. / 5.)*m*R*R;
}
void Cuerpo::BorreFuerza_Torque(void) {
	F.cargue(0, 0, 0);  tau.cargue(0, 0, 0);
}
void Cuerpo::AgregueFuerza(vector3D F0) {
	F += F0;
}
void Cuerpo::AgregueTorque(vector3D tau0) {
	tau += tau0;
}
void Cuerpo::Mueva_r(double dt, double Constante) {
	r += V*(Constante*dt);  theta += omega.z()*Constante*dt;
}
void Cuerpo::Mueva_V(double dt, double Constante) {
	V += F*(Constante*dt / m); omega += tau*(Constante*dt / I);
}
void Cuerpo::Dibujese(void) {
	cout << ", " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t), "
		<< r.x() << "+" << R*cos(theta) / 7.0 << "*t," << r.y() << "+" << R*sin(theta) / 7.0 << "*t";
}
//----------------------Clase Colisionador---------------
class Colisionador {
private:
	vector3D ele[N + 7][N + 7]; bool EstoyEnColision[N + 7][N + 7];
public:
	void Inicie(void);
	void CalculeTodasLasFuerzas(Cuerpo* Grano, double dt, int Shake);
	void CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,
		vector3D & ele, bool & EstoyEnColision, double dt);
};
void Colisionador::Inicie(void) {
	int i, j;
	for (i = 0;i<N + 1;i++)
		for (j = i + 1;j<N + 7;j++) {
			ele[i][j].cargue(0, 0, 0); EstoyEnColision[i][j] = false;
		}
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano, double dt, int Shake) {
	int i, j;
	vector3D g_vector; g_vector.cargue(0, 0, -g);
	if (boundary == 0) {
		//Borrar todas las fuerzas y torques
		for (i = 0;i < N + 7;i++) Grano[i].BorreFuerza_Torque();
		//Calcular todas las fuerzas entre parejas de planetas
		for (i = 0;i < N + 1;i++)
			for (j = i + 1;j < N + 7;j++)
				CalculeLaFuerzaEntre(Grano[i], Grano[j], ele[i][j], EstoyEnColision[i][j], dt);
		for (i = 0;i < N + 1;i++) { Grano[i].AgregueFuerza(g_vector*Grano[i].m); }
		//Borrar fuerzas sobre las paredes
		for (i = N + 1;i < N + 7;i++) Grano[i].BorreFuerza_Torque();
		//Fuerza armonica para agitar la caja
		vector3D z_1;
		z_1.cargue(0, 0, 1);
		if (Shake == 1) {
			Grano[N + 1].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
			Grano[N + 2].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
			Grano[N + 3].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
			Grano[N + 4].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
			Grano[N + 5].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
			Grano[N + 6].AgregueFuerza(-k_Shake*(Grano[N + 5].Getz() + (Lz / 2 + Grano[N + 5].R))*z_1);
		}
	}
	else {
		//Borrar todas las fuerzas y torques
		for (i = 0;i<N + 3;i++) Grano[i].BorreFuerza_Torque();
		//Calcular todas las fuerzas entre parejas de planetas
		for (i = 0;i<N + 1;i++)
			for (j = i + 1;j<N + 3;j++)
				CalculeLaFuerzaEntre(Grano[i], Grano[j], ele[i][j], EstoyEnColision[i][j], dt);
		for (i = 0;i<N + 1;i++) { Grano[i].AgregueFuerza(g_vector*Grano[i].m); }
		//Borrar fuerzas sobre las paredes
		for (i = N + 1;i < N + 3;i++) Grano[i].BorreFuerza_Torque();
		//Fuerza armonica para agitar la caja
		vector3D z_1;
		z_1.cargue(0, 0, 1);
		if (Shake == 1) {
			Grano[N + 1].AgregueFuerza(-k_Shake*(Grano[N + 1].Getz() + (Lz / 2 + Grano[N + 1].R))*z_1);
			Grano[N + 2].AgregueFuerza(-k_Shake*(Grano[N + 1].Getz() + (Lz / 2 + Grano[N + 1].R))*z_1);
		}
		//Cilindrical boundary
		vector3D r_normal, r_tangencial, F, Fn, Ft, Vc, Vcn, Vct;
		double componenteFn, componenteVcn, normaVct, normaFt, Ftmax;
		double xp, yp, s;
		double ERFF = 1e-8;

		for (i = 0;i < N + 1;i++) {
			xp = Grano[i].Getx();  yp = Grano[i].Gety();  s = (sqrt(xp*xp + yp*yp) + Grano[i].R) - R_cil;
			if (s > 0) {
				r_normal.cargue(-xp / sqrt(xp*xp + yp*yp), -yp / sqrt(xp*xp + yp*yp), 0);
				//Calcular velocidad de contacto y vector tangente
				Vc = (Grano[i].V) - (Grano[i].omega ^ (r_normal*Grano[i].R));
				componenteVcn = Vc*r_normal; Vcn = r_normal*componenteVcn; Vct = Vc - Vcn; normaVct = norma(Vct);

				//FUERZAS NORMALES
				//Fuerza Hertz
				componenteFn = K*pow(s, 1.5);
				//Dispacion plastica
				componenteFn -= Grano[i].m*sqrt(s)*Gamma*componenteVcn; if (componenteFn < 0) componenteFn = 0;
				Fn = r_normal*componenteFn;

				//FUERZAS TANGENCIALES
				//Fuerza estatica
				r_tangencial += (Vct*dt);
				Ft = r_tangencial*(-Kcundall);
				//Fuerza cinetica
				Ftmax = MU*componenteFn;  normaFt = norma(Ft);
				if (normaFt > Ftmax) Ft = r_tangencial*(-Ftmax / norma(r_tangencial));
				//Construir fuerza total
				F = Fn + Ft;
				Grano[i].AgregueFuerza(F);      Grano[i].AgregueTorque((r_normal*(-Grano[i].R)) ^ Ft);
			}
		}
	}
}
void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,
	vector3D & ele, bool & EstoyEnColision, double dt) {
	vector3D F2, Fn, Ft;
	vector3D Vc, Vcn, Vct;
	vector3D r21, n, t;

	double R1, R2, d21, s;
	double m1, m2, m21;
	double componenteVcn, normaVct, componenteFn;
	double normaFt, Ftmax;

	double ERFF = 1e-8;
	r21 = Grano2.r - Grano1.r; d21 = norma(r21); s = (Grano1.R + Grano2.R) - d21;
	if (s>0) { //Si se chocan,
			   //Geomatria y dinamica del contacto
		m1 = Grano1.m;   R1 = Grano1.R;
		m2 = Grano2.m;   R2 = Grano2.R;
		m21 = (m1*m2) / (m1 + m2);
		n = r21 / d21;
		//Calcular velocidad de contacto y vector tangente
		Vc = (Grano2.V - Grano1.V) - (Grano2.omega ^ (n*R2)) - (Grano1.omega ^ (n*R1));
		componenteVcn = Vc*n; Vcn = n*componenteVcn; Vct = Vc - Vcn; normaVct = norma(Vct);
		if (normaVct<ERFF) t.cargue(0, 0, 0); else t = Vct / normaVct;

		//FUERZAS NORMALES
		//Fuerza Hertz
		componenteFn = K*pow(s, 1.5);
		//Dispacion plastica
		componenteFn -= m21*sqrt(s)*Gamma*componenteVcn; if (componenteFn<0) componenteFn = 0;
		Fn = n*componenteFn;

		//FUERZAS TANGENCIALES
		//Fuerza estatica
		ele += (Vct*dt);
		Ft = ele*(-Kcundall);
		//Fuerza cinetica
		Ftmax = MU*componenteFn;  normaFt = norma(Ft);
		if (normaFt>Ftmax) Ft = ele*(-Ftmax / norma(ele));

		//Construir fuerza total
		F2 = Fn + Ft;
		Grano2.AgregueFuerza(F2);      Grano2.AgregueTorque((n*(-R2)) ^ Ft);
		Grano1.AgregueFuerza(F2*(-1)); Grano1.AgregueTorque((n*R1) ^ (Ft*(-1)));

		EstoyEnColision = true;
	}
	else if (EstoyEnColision == true) {
		ele.cargue(0, 0, 0); EstoyEnColision = false;
	}
}
//----------------------Funciones Globales---------------
void Init_POV(Cuerpo * Grano, double Number_of_vectors) {
	ofstream Granular_pov;
	ofstream Granular_ini;
	Granular_pov.open("./Output/Granular.pov");
	Granular_pov << "#include \"colors.inc\"\n";

	Granular_pov << "#declare cubical_box=\n";
	Granular_pov << " 	union{\n";
	Granular_pov << " 	box{ <" << -Lx / 2 << "," << -Ly / 2 << "," << -(1. + Amplitud)*Lz / 2 << ">,<" << Lx / 2 << "," << -Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << "> texture{pigment{ color rgbt<1,1,1,0.8>}  finish { phong 1.0 }} }   \n";
	Granular_pov << " 	box{ <-" << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<-" << Lx / 2 << "," << Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << "> texture{pigment{ color rgbt<1,1,1,0.8>}  finish { phong 1.0 }} }  \n";
	Granular_pov << " 	box{ <-" << Lx / 2 << "," << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<" << Lx / 2 << "," << Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << "> texture{pigment{ color rgbt<1,1,1,0.8>}  finish { phong 1.0 }} }  \n";
	Granular_pov << " 	box{ <" << Lx / 2 << "," << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<" << Lx / 2 << ",-" << Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << "> texture{pigment{ color rgbt<1,1,1,0.8>}  finish { phong 1.0 }} } \n";
	Granular_pov << " 	box{ <-" << Lx / 2 << ",-" << Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << ">,<" << Lx / 2 << "," << Ly / 2 << "," << (1. + Amplitud)*Lz / 2 << "> texture{pigment{ color rgbt<1,1,1,0.8>}  finish { phong 1.0 }} } \n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<-" << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<-" << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ < " << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ < " << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";

	Granular_pov << "	cylinder{ <-" << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,<-" << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,<-" << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ < " << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ < " << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";

	Granular_pov << "	cylinder{ <-" << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ",-" << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ",-" << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ", " << Ly / 2 << ",-" << (1. + Amplitud)*Lz / 2 << ">,0.3}\n";
	Granular_pov << "	cylinder{ <-" << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,< " << Lx / 2 << ", " << Ly / 2 << ", " << (1. + Amplitud)*Lz / 2 << ">,0.3}}\n";

	Granular_pov << "#declare cilindrical_box=\n";
	Granular_pov << "cylinder{ <0,0,-" << (1. + Amplitud)*Lz / 2 << ">,<0,0," << (1.)*Lz / 2 << ">," << R_cil << "\n";
	Granular_pov << "        texture{ pigment{ color rgbt<0,1,1,0.8>}\n";
	Granular_pov << "                finish { phong 1.0}\n";
	Granular_pov << "                }\n";
	Granular_pov << "      rotate<0,0,0> translate<0,0,0>\n";
	Granular_pov << "    }\n";

	Granular_pov << "global_settings {assumed_gamma 1.0} \n";
	Granular_pov << "camera{";
	Granular_pov << "      right x*image_width/image_height \n";
	Granular_pov << "      location  <0.0 , " << 2.*Lz << " ,0.0>\n";
	Granular_pov << "      look_at   <0.0 , 0.0 , 0.0> }light_source{ <1500,2500,-2500> \n";
	Granular_pov << "      color rgb<1,1,1> } \n";
	Granular_pov << "      sky_sphere{ pigment{color rgb<1,1,1>}} \n";
	Granular_pov << "#declare Vector = array[" << Number_of_vectors << "];\n";
	Granular_pov << "#fopen MyFile  \"data.txt\"   read\n";
	Granular_pov << "#local i = 0;\n";
	Granular_pov << "#while(i<" << Number_of_vectors << ")\n";
	Granular_pov << "    #local Vector[i]=<0,0,0>;\n";
	Granular_pov << "    #read (MyFile,Vector[i])\n";
	Granular_pov << "    #local i = i + 1;\n";
	Granular_pov << "    #end\n";
	Granular_pov << "#declare Grano=\n";
	Granular_pov << "    sphere{<0,0,0>, " << Grano[0].Get_R() << "\n";
	Granular_pov << "    material{texture { pigment{ color rgbt<0,1,1,0>} \n";
	Granular_pov << "            finish { phong 1.0 }}} }\n";
	Granular_pov << "#declare Nut=\n";
	Granular_pov << "    sphere{<0,0,0>, " << Grano[N].Get_R() << "\n";
	Granular_pov << "    material{texture { pigment{ color rgbt<1,1,1,0>} \n";
	Granular_pov << "            finish { phong 1.0 }}} }\n";
	for (int i = 0;i<N;i++) {
		Granular_pov << "object{Grano  translate Vector[" << i << "+(frame_number-1)*" << N + 2 << "]}\n";
	}
	Granular_pov << "object{Nut  translate Vector[" << N << "+(frame_number-1)*" << N + 2 << "]}\n";
	if (boundary == 0) {
		Granular_pov << "object{cubical_box}\n";
	}
	else {
		Granular_pov << "object{cilindrical_box translate Vector[" << N + 1 << "+(frame_number-1)*" << N + 2 << "]}\n";
	}
	Granular_pov.close();
	//----------------------------------------------------------
	Granular_ini.open("./Output/Granular.ini");
	Granular_ini << "; POV-Ray animation ini file\n";
	Granular_ini << "Output_File_Type=J\n";
	Granular_ini << "Frame_Step = 500\n";
	Granular_ini << "Antialias=Off\n";
	Granular_ini << "Antialias_Threshold=0.1\n";
	Granular_ini << "Antialias_Depth=2\n";
	Granular_ini << "Input_File_Name=\"Granular.pov\"\n";
	Granular_ini << "Initial_Frame=1\n";
	Granular_ini << "Final_Frame=" << Number_of_vectors / (N + 2) << "\n";
	Granular_ini << "Initial_Clock=0\n";
	Granular_ini << "Final_Clock=1\n";
	Granular_ini << "Cyclic_Animation=on\n";
	Granular_ini << "Pause_when_Done=off\n";
	Granular_ini.close();

}
//--------------------Programa Principal------------------
int main(void) {
	double t, dt = 1e-3;
	Cuerpo Grano[N + 7]; int i, j, k;
	Colisionador Newton;
	Crandom ran64(1); double theta, phi;

	double m0 = 2.5, R0 = 1.8, V = 100., omega0 = 0;
	double Rpared = 10000, Mpared = 1000;

	double R_Nut = 6.3, M_Nut = m0;

	double T = Lx / V, t_Shake = 20*T, tmax = 200*T;

	double dx = Lx / (Nx + 1), dy = Ly / (Ny + 1), dz = 2*R0;
	double dr = R_cil / (Nr + 1), dtheta = 2 * M_PI / (Ntheta);

	ofstream data_pov;
	ofstream data_z;
	data_pov.open("./Output/data.txt");
	data_z.open("./Output/z_data.txt");
	int Shake = 0;
	if (boundary == 0) {
		//Init all particles
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {
					phi = 2 * M_PI*ran64.r();
					theta = M_PI*ran64.r();
					Grano[i + Ny*j + Nx*Ny*k].Inicie(-(Lx / 2) + (i + 1)*dx, -(Ly / 2) + (j + 1)*dy, (-0.5*Lz*(1+Amplitud)+2*R_Nut)+(k + 1)*dz, 0, 0, 0, 0, 0, 0, omega0, m0, R0);
				}
			}
		}
		//Nut
		Grano[N].Inicie(0, 0, -(1. + Amplitud)*(0.5*Lz) + R_Nut, 0, 0, 0, 0, 0, 0, 0, M_Nut, R_Nut);
		//Pared arriba
		Grano[N + 1].Inicie(0, (Ly / 2) + Rpared, 0, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
		//Pared abajo
		Grano[N + 2].Inicie(0, -(Ly / 2) - Rpared, 0, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
		//Pared derecha
		Grano[N + 3].Inicie((Lx / 2) + Rpared, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
		//Pared izquierda
		Grano[N + 4].Inicie(-(Lx / 2) - Rpared, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);

		Grano[N + 5].Inicie(0, 0, -(1. + Amplitud)*(Lz / 2) - Rpared, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);

		Grano[N + 6].Inicie(0, 0, (Lz / 2) + Rpared, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);

		Init_POV(Grano, (N + 2)*tmax / dt);
		//Muevase con Omelian
		for (t = 0;t<tmax;t += dt) {
			//Generate data
			//Pariculas Pov
			for (int i = 0;i<N + 2;i++) { data_pov << "<" << (Grano[i].Getx()) << "," << (Grano[i].Gety()) << "," << (Grano[i].Getz()) << ">" << ","; }
			//Altura Nuez
			data_z << t << " " << Grano[N].Getz() << endl;
			if (t>t_Shake) { Shake = 1; };
			//Muevase con Omelyan PEFRL
			for (i = 0;i<N + 7;i++) Grano[i].Mueva_r(dt, Zeta);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 7;i++) Grano[i].Mueva_V(dt, (1 - 2 * Lambda) / 2);
			for (i = 0;i<N + 7;i++) Grano[i].Mueva_r(dt, Xi);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 7;i++) Grano[i].Mueva_V(dt, Lambda);
			for (i = 0;i<N + 7;i++) Grano[i].Mueva_r(dt, 1 - 2 * (Xi + Zeta));
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 7;i++) Grano[i].Mueva_V(dt, Lambda);
			for (i = 0;i<N + 7;i++) Grano[i].Mueva_r(dt, Xi);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 7;i++) Grano[i].Mueva_V(dt, (1 - 2 * Lambda) / 2);
			for (i = 0;i<N + 7;i++) Grano[i].Mueva_r(dt, Zeta);
		}
	}
	else {
		//Init all particles
		for (k = 0; k < Nz; k++) {
			for (j = 0; j < Ntheta; j++) {
				for (i = 0; i < Nr; i++) {
					Grano[i + j*Nr + Nr*Ntheta*k].Inicie((i + 1)*dr*cos(dtheta*j), (i + 1)*dr*sin(dtheta*j), (-0.5*Lz*(1+Amplitud)+2*R_Nut)+(k + 1)*dz, 0, 0 , 0, 0, 0, 0, omega0, m0, R0);
				}
			}
		}
		//Nut
		Grano[N].Inicie(0, 0, -(1.+ Amplitud)*(0.5*Lz) + R_Nut, 0, 0, 0, 0, 0, 0, 0, M_Nut, R_Nut);
		//Pared arriba
		Grano[N + 1].Inicie(0, 0, -(1. + Amplitud)*(Lz / 2) - Rpared, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
		//Pared abajo
		Grano[N + 2].Inicie(0, 0, (Lz / 2) + Rpared, 0, 0, 0, 0, 0, 0, 0, Mpared, Rpared);

		Init_POV(Grano, (N + 2)*tmax / dt);
		//Muevase con Omelian
		for (t = 0;t<tmax;t += dt) {
			//Generate data
			//Pariculas Pov
			for (int i = 0;i<N + 1;i++) { data_pov << "<" << (Grano[i].Getx()) << "," << (Grano[i].Gety()) << "," << (Grano[i].Getz()) << ">" << ","; }
			data_pov << "<" << (Grano[N + 1].Getx()) << "," << (Grano[N + 1].Gety()) << "," << (Grano[N + 1].Getz() + 0.5*Lz*(1.+ Amplitud) + Rpared) << ">" << ",";
			//Altura Nuez
			data_z << t << " " << Grano[N].Getz() << endl;
			if (t>t_Shake){ Shake = 1; };
			//Muevase con Omelyan PEFRL
			for (i = 0;i<N + 3;i++) Grano[i].Mueva_r(dt, Zeta);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 3;i++) Grano[i].Mueva_V(dt, (1 - 2 * Lambda) / 2);
			for (i = 0;i<N + 3;i++) Grano[i].Mueva_r(dt, Xi);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 3;i++) Grano[i].Mueva_V(dt, Lambda);
			for (i = 0;i<N + 3;i++) Grano[i].Mueva_r(dt, 1 - 2 * (Xi + Zeta));
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 3;i++) Grano[i].Mueva_V(dt, Lambda);
			for (i = 0;i<N + 3;i++) Grano[i].Mueva_r(dt, Xi);
			Newton.CalculeTodasLasFuerzas(Grano, dt, Shake); for (i = 0;i<N + 3;i++) Grano[i].Mueva_V(dt, (1 - 2 * Lambda) / 2);
			for (i = 0;i<N + 3;i++) Grano[i].Mueva_r(dt, Zeta);
		}
	}
	data_pov.close();
	data_z.close();
	return 0;
}