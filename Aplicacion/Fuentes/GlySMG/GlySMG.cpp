/*
	numero de caso
	mono/multi tarea (0 o 1
	intervalo de tiempo en años
	geometria (E o C)
	lado del cubo
	generar velocidades (velocidad en al/s o H para generarlas con la ley de Hubble, H10 = módulo x 10)
	distancia de colision en al
	radio de corte en 'al'
	frecuencia para renovar la tabla de vecinas (0 para no usar tabla)
	radio de la tabla de vecinas en 'al'
	datos de la muestra (senda y nombre de fichero, R (en retícula) o A (azar)
	pasos a simular
	pasos para monitorizar
	numero de hilos
	numero de particulas
	hasta un máximo de 32, puede no haber ninguna
	numero de particula a seguir (de 1 hasta numero de particulas)
	numero de particula a seguir (de 1 hasta numero de particulas)
	...

 0 1 0.000000031687536 E 1 H 1E-13 2 0 1 F:/Gly/solar.csv 31558150 86400 1 1 3 4 5 11
 0 0 0.000000031687536 E 1 H 1E-13 2 0 1 F:/Gly/solar.csv 31558150 86400 1 1 1 2 3 4 5 6 7 8 9 10 11 12
 101 1 0.00000003169 E 1 H 1E-13 2 0 1 F:/Gly/smg_sal_AVX256_L_m1_1h_E_12p_C100_1.csv 31557600 3155760 1 1
 100 1 1000 E 180000000 H 5E+04 90000000 100 100000000 R 4000000 10000 30 10680

 Densidad de referencia 2,1E-22 particulas / al^3
 21.960 partículas -> lado cubo = 4,716E08, radio esfera = 2,926E08 diametro = 5,581E08

 0 1 1000 C 4.716E08 H 5E+04 1.5E08 100 2.0E08 R 1E06 1E04 30 21960 1 7320 14640 21952

 0 1 10000 E 5.581E08 1.2E-13 5E04 1.9E08 100 2.3E08 A 1E06 1E03 30 21960 1 7320 14640 21960

 0 1 100 E 5E06 H1E03 5E01 2.1E08 0 2.3E08 A 1E06 1E03 30 21960

 Una simulación de 4.0E09 años expande el radio aproximadamente en 50.0E06 'al', por lo tanto partimos de él reducido

 0 1 1000 E 4.581E08 H 5E+04 1.5E08 100 2.0E08 R 1E06 1E04 30 21960 1 7320 14640 20672
 0 1 1000 C 3.716E08 H 5E+04 1.5E08 100 2.0E08 R 1E06 1E04 30 21960 1 7320 14640 21952

 */

#include "Header.h"

using namespace std;
using namespace chrono;

void CapacidadesCPU();
void EsperaCaracter();
bool BuscaNombreDeProceso(DWORD processID);
void SolicitaCaso();
void SolicitaTipoProceso();
void SolicitaPasoT();
void SolicitaGeometria();
void SolicitaLadoCubo();
void SolicitaVelocidadInicial();
void SolicitaDisColision();
void SolicitaRadioCorte();
void SolicitaFrecuenciaVecinas();
void SolicitaRadioVecinas();
void SolicitaGalaxias();
char* SolicitaFichero(char* s);
void SolicitaDuracion();
void SolicitaFrecuenciaMonitor();
void SolicitaHilos();
bool CompletaParametros();
char* SeparaMiles(char* s, unsigned long n);
char* ConvPtoComa(char* s);
int DatoColumna(char* linea, int i, char* campo);
int LeeInicio(char* fichero);
double Covarianza(int n, double x[], double y[]);
bool GeneraMasasAzar(double vmin, double vmax);
bool GeneraCoordenadasAzar(double lcubo, bool esfera);
bool GeneraCoordenadasReticula(double lcubo, bool esfera);
bool GeneraVelocidadesAzar(bool radial, double vmax, double lcubo);
bool GeneraVelocidadesHubble(double lcubo, double vel_inicial);
bool CalculaVelocidadesIntermediasIniciales();
bool Arranque(char* ficheromuestra, unsigned long hasta);
void Seguimiento(unsigned long iteracion);
double MediaVecinas();
void MuestraDatosParticula(FILE* ef, int i);
void ImprimeDatosParticula(FILE* ef, int i);
void ImprimeSeguimiento(unsigned long n);
void Terminacion(unsigned long iteracion, unsigned long hasta, const char* rotulo);
void ProcesoMono(char* ficheromuestra, unsigned long hasta);
void PreparaMP(char* ficheromuestra, unsigned long hasta);
unsigned long ProcesoMP(unsigned long hasta, int gh);
void Aceleraciones();
void AceleracionesTramo(int hilo, int desde, int hasta);
bool IncrementaTMono();
bool EscribeDatos(int n);
void CentroMasas(double* cm);
double DistanciaMediaAcm(double* cm);
double DesviacionEstandarDisAcm(double* cm, double dm);
double Ecinetica();
double EcineticaParticula(int i);
double Epotencial(double* epr, double* ndentro);
double EpotencialParticula(int i, double* epr);
void TrataColision(unsigned long iteracion, int i, int j, double* pd2, double* pd);
void TrataColisionHilo(unsigned long iteracion, int i, int j, double* pd2, double* pd);
void ActualizaTablaVecinas();
unsigned int WINAPI HiloAceleraciones(LPVOID param);
unsigned int WINAPI HiloVecinas(LPVOID param);
void LiberaMemoria();
bool AsignaMemoria(int nparticulas);
//bool FusionaMono(unsigned long iteracion, int i, int j);
//bool Fusiona(unsigned long iteracion, int i, int j);

BOOL CtrlHandler(DWORD fdwCtrlType);

#ifdef UNIDADES_SMI
// Se trabaja en unidades del sistema métrico internacional, aunque los resultados se muestran en 'al'

const double AL = 9460730472580800;						// m    (9,461e+15 m)
const double G = 6.67408E-11;							// N m^2 / kg^2  =  m^3/ (kg s^2)
const double MASA_VIA_LACTEA = 2.0E+42;					// kg 
#else
const double AL = 1.0;
const double G = 1.567267014756E-28;
const double MASA_VIA_LACTEA = 1.00578325370883E+12;
#endif // UNIDADES_SMI

const double H0 = 2.26963E-18;							// 70 (km/s)/Mpc = 2,26963E-18 (m/s)/m
const double SEGUNDOS_ANNO = 365.25636 * 24.0 * 3600.0;

const int MAX_PARTICULAS = 22000;						// el limite depende del producto MAX_PARTICULAS * 2 * MAX_VECINAS
const int MAX_VECINAS = 22000;
const int MAX_HILOS = 64;
const int NUM_HISTO = 30;
const int MAX_SEGUIR = 32;
const char filog[] = "F:\\Gly\\log\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_log.csv";
const char ficol[] = "F:\\Gly\\log\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_col.csv";
const char fisal[] = "F:\\Gly\\sal\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_sal_%d.csv";
const char fiseg[] = "F:\\Gly\\seg\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_seg_%d.csv";
const char fihtc[] = "F:\\Gly\\htg\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_htc.csv";
const char fihtv[] = "F:\\Gly\\htg\\smg___%s%s_m%d_%dh_%s_%s_%dp_caso%d%s_htv.csv";
char rot_md[3];
char rot_avx[6];
char rot_su[5];

int n_hilos;
HANDLE hilos[MAX_HILOS];
HANDLE semaforo_hilo[MAX_HILOS];
int hilos_terminados;
HANDLE semaforo_itera;
CRITICAL_SECTION seccion_fin_hilos;
CRITICAL_SECTION seccion_colision;
#ifdef HILO_VECINAS
HANDLE hilo_vecinas;
HANDLE semaforo_vecinas;
#endif // HILO_VECINAS

int n_caso;
char tipo_inicio[4];
int t_proceso;
double paso_t;
double medio_paso_t;
double Gxpaso_t;
double Gxmedio_paso_t;
bool esfera;
char geometria[2];
double lado_cubo;
int modo_v;
double vel_inicial;
double dis_colision;
double dis_critica;
double radio_corte;
double radio_corte2;
unsigned long f_vecinas;
int u_vecinas;
int t_vecinas_enuso;
int t_vecinas_actualizada;
double radio_vecinas;
double radio_vecinas2;
int max_absoluto_vecinas;
unsigned long max_iteraciones;
unsigned long iteracion;
unsigned long frecu_monitorizar;
unsigned long frecu_salidas;
int nparticulas;
int nsalida;
int ncolisiones;
int ncriticas;
int nficticias;
int nseguir;
int particulas_seguir[MAX_SEGUIR];
double coor_ini[3][MAX_SEGUIR];
double recorrido[MAX_SEGUIR];
FILE* eseg[MAX_SEGUIR];
#ifdef AVX2
const int bloque_mt = 4;
int npt_mt;
int npm4_mt;
int nmj_mt;
__m256d pasot_mt;
__m256d gxpasot_mt;
__m256d gxmediopasot_mt;
__m256d gxmediopasot_menos_mt;
__m256d coor_antes_mt;
__m256d coor_despues_mt;
__m256d vi_antes_mt;
__m256d vi_despues_mt;
__m256d v_despues_mt;
__m256d a_antes_mt;
double* val_coor_mt;
double* val_vi_mt;
double* val_v_mt;
#endif // AVX2

FILE* elog;
FILE* ecol;

char fto1[128];
char fto2[128];

#ifdef MEMORIA_DINAMICA
double* masas = NULL;
double** coordenadas = NULL;
double** velocidades = NULL;
double** velintermedias = NULL;
double** aceleraciones = NULL;
int** n_vecinas = NULL;
#else
// alignas(32) !! Aumenta el tiempo de cálculo !!

double masas[MAX_PARTICULAS];
double coordenadas[3][MAX_PARTICULAS];
double velocidades[3][MAX_PARTICULAS];
double velintermedias[3][MAX_PARTICULAS];
double aceleraciones[3][MAX_PARTICULAS];
int n_vecinas[MAX_PARTICULAS][2];
//int vecinas[MAX_PARTICULAS][2][MAX_VECINAS];
#endif // MEMORIA_DINAMICA
int*** vecinas = NULL;

time_t reloji;
time_t reloja;
double segundos;
bool CANCELADO;

bool HW_MMX;
bool HW_x64;
bool HW_ABM;      // Advanced Bit Manipulation
bool HW_RDRAND;
bool HW_BMI1;
bool HW_BMI2;
bool HW_ADX;
bool HW_PREFETCHWT1;

//  SIMD: 128-bit

bool HW_SSE;
bool HW_SSE2;
bool HW_SSE3;
bool HW_SSSE3;
bool HW_SSE41;
bool HW_SSE42;
bool HW_SSE4a;
bool HW_AES;
bool HW_SHA;

//  SIMD: 256-bit

bool HW_AVX;
bool HW_XOP;
bool HW_FMA3;
bool HW_FMA4;
bool HW_AVX2;

//  SIMD: 512-bit

bool HW_AVX512F;    //  AVX512 Foundation
bool HW_AVX512CD;   //  AVX512 Conflict Detection
bool HW_AVX512PF;   //  AVX512 Prefetch
bool HW_AVX512ER;   //  AVX512 Exponential + Reciprocal
bool HW_AVX512VL;   //  AVX512 Vector Length Extensions
bool HW_AVX512BW;   //  AVX512 Byte + Word
bool HW_AVX512DQ;   //  AVX512 Doubleword + Quadword
bool HW_AVX512IFMA; //  AVX512 Integer 52-bit Fused Multiply-Add
bool HW_AVX512VBMI; //  AVX512 Vector Byte Manipulation Instructions

void CapacidadesCPU() {

	int info[4];
	cpuid(info, 0);
	int nIds = info[0];
	cpuid(info, 0x80000000);
	unsigned nExIds = info[0];

	//  Detect Features

	if (nIds >= 0x00000001) {
		cpuid(info, 0x00000001);
		HW_MMX = (info[3] & ((int)1 << 23)) != 0;
		HW_SSE = (info[3] & ((int)1 << 25)) != 0;
		HW_SSE2 = (info[3] & ((int)1 << 26)) != 0;
		HW_SSE3 = (info[2] & ((int)1 << 0)) != 0;

		HW_SSSE3 = (info[2] & ((int)1 << 9)) != 0;
		HW_SSE41 = (info[2] & ((int)1 << 19)) != 0;
		HW_SSE42 = (info[2] & ((int)1 << 20)) != 0;
		HW_AES = (info[2] & ((int)1 << 25)) != 0;

		HW_AVX = (info[2] & ((int)1 << 28)) != 0;
		HW_FMA3 = (info[2] & ((int)1 << 12)) != 0;

		HW_RDRAND = (info[2] & ((int)1 << 30)) != 0;
	}
	if (nIds >= 0x00000007) {
		cpuid(info, 0x00000007);
		HW_AVX2 = (info[1] & ((int)1 << 5)) != 0;

		HW_BMI1 = (info[1] & ((int)1 << 3)) != 0;
		HW_BMI2 = (info[1] & ((int)1 << 8)) != 0;
		HW_ADX = (info[1] & ((int)1 << 19)) != 0;
		HW_SHA = (info[1] & ((int)1 << 29)) != 0;
		HW_PREFETCHWT1 = (info[2] & ((int)1 << 0)) != 0;

		HW_AVX512F = (info[1] & ((int)1 << 16)) != 0;
		HW_AVX512CD = (info[1] & ((int)1 << 28)) != 0;
		HW_AVX512PF = (info[1] & ((int)1 << 26)) != 0;
		HW_AVX512ER = (info[1] & ((int)1 << 27)) != 0;
		HW_AVX512VL = (info[1] & ((int)1 << 31)) != 0;
		HW_AVX512BW = (info[1] & ((int)1 << 30)) != 0;
		HW_AVX512DQ = (info[1] & ((int)1 << 17)) != 0;
		HW_AVX512IFMA = (info[1] & ((int)1 << 21)) != 0;
		HW_AVX512VBMI = (info[2] & ((int)1 << 1)) != 0;
	}
	if (nExIds >= 0x80000001) {
		cpuid(info, 0x80000001);
		HW_x64 = (info[3] & ((int)1 << 29)) != 0;
		HW_ABM = (info[2] & ((int)1 << 5)) != 0;
		HW_SSE4a = (info[2] & ((int)1 << 6)) != 0;
		HW_FMA4 = (info[2] & ((int)1 << 16)) != 0;
		HW_XOP = (info[2] & ((int)1 << 11)) != 0;
	}
}

int subversion = 15;
int version = 4;

int main(int argc, char* argv[])
{
	CapacidadesCPU();
#ifdef AVX2
	if (!HW_AVX2) {
		printf("Esta versión de la aplicación necesita AVX2 y la CPU no la admite\n");
		EsperaCaracter();
		return -1;
	}
	if (HW_AVX512F) {
		printf("Aviso. Esta CPU admite AVX512 y la aplicacion no la utiliza\n");
	}
#else
	if (HW_AVX512F) {
		printf("Aviso. Esta CPU admite AVX512 y la aplicacion no la utiliza\n");
	}
	else {
		if (HW_AVX2) {
			printf("Aviso. Esta CPU admite AVX2 y la aplicacion no la utiliza\n");
		}
	}
#endif
	// Para que en los números la coma sea la coma decimal y el punto el punto de millares (localización española)

	locale loc_nueva("");
	locale loc_antes = locale::global(loc_nueva);
	/*cout << "El nombre de la localizacion anterior era: " << loc_antes.name() << "." << endl;
	cout << "El nombre de la localizacion actual es   : " << loc_nueva.name() << "." << endl;*/

	HWND console = GetConsoleWindow();
	RECT r = { 2, 2, 1880, 600 };
	MoveWindow(console, r.left, r.top, r.right, r.bottom, TRUE);
	char titulo[128];
	TCHAR ctitulo[128];
	sprintf_s(titulo, 128, "Simulacion Movimiento de particulas Gravitatorias. v %d.%02d", version, subversion);
	_stprintf_s(ctitulo, sizeof(ctitulo) / sizeof(TCHAR), _T("Simulacion Movimiento de particulas Gravitatorias. v %d.%02d"), version, subversion);
	SetConsoleTitle(ctitulo);

	if (!SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE))
	{
		printf("\nERROR [%d] al establecer el controlador de eventos\n", GetLastError());
		EsperaCaracter();
		return -1;
	}
	int resultado;
	char carpeta[17];
	sprintf_s(carpeta, 16, "F:\\Gly\\log");
	resultado = _mkdir(carpeta);
	sprintf_s(carpeta, 16, "F:\\Gly\\sal");
	resultado = _mkdir(carpeta);
	sprintf_s(carpeta, 16, "F:\\Gly\\seg");
	resultado = _mkdir(carpeta);
	sprintf_s(carpeta, 16, "F:\\Gly\\htg");
	resultado = _mkdir(carpeta);
	t_vecinas_actualizada = t_vecinas_enuso = -1;
	max_absoluto_vecinas = 0;
	char ficheromuestra[512] = { 0 };
	if (argc > 15) {
		n_caso = atoi(argv[1]);
		t_proceso = atoi(argv[2]);
		paso_t = atof(ConvPtoComa(argv[3]));
		esfera = argv[4][0] == 'E';
		lado_cubo = atof(ConvPtoComa(argv[5]));
		if (argv[6][0] == 'H' || argv[6][0] == 'h') {

			// Módulo Hubble y dirección radial 

			modo_v = 0;
			if (strlen(argv[6]) == 1) {
				vel_inicial = 1.0;
			}
			else {
				argv[6][0] = '0';
				vel_inicial = atof(ConvPtoComa(argv[6]));
			}
		}
		else {

			// Módulo al azar

			vel_inicial = atof(ConvPtoComa(argv[6]));
			if (vel_inicial >= 0.0) {

				// Dirección al azar

				modo_v = 1;
			}
			else {

				// Dirección radial

				modo_v = 2;
			}
		}
		dis_colision = atof(ConvPtoComa(argv[7]));
		radio_corte = atof(ConvPtoComa(argv[8]));
		f_vecinas = atol(argv[9]);
		radio_vecinas = atof(ConvPtoComa(argv[10]));
		strcpy_s(ficheromuestra, 512, argv[11]);
		max_iteraciones = (unsigned long)atof(argv[12]);
		frecu_monitorizar = (unsigned long)atof(argv[13]);
		n_hilos = atoi(argv[14]);
		if (n_hilos > MAX_HILOS) {
			printf("\nDemasiados hilos (maximo = %d)\n", MAX_HILOS);
			EsperaCaracter();
			return -1;
		}
		nparticulas = atoi(argv[15]);
		if (nparticulas > MAX_PARTICULAS) {
			printf("\nDemasiadas partículas (maximo = %d)\n", MAX_PARTICULAS);
			EsperaCaracter();
			return -1;
		}
		if (argc > 16) {
			nseguir = argc - 16;
			if (nseguir > MAX_SEGUIR)nseguir = MAX_SEGUIR;
			for (int i = 0; i < nseguir; i++) {
				particulas_seguir[i] = atoi(argv[16 + i]) - 1;
			}
		}
		else {
			nseguir = 1;
			particulas_seguir[0] = 0;
		}
		if (!CompletaParametros()) {
			EsperaCaracter();
			return -1;
		}
		if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'A' || ficheromuestra[0] == 'a')) {
			if (modo_v == 0) {
				sprintf_s(tipo_inicio, 4, "%s", "AH");
			}
			else if (modo_v == 1) {
				sprintf_s(tipo_inicio, 4, "%s", "AA");
			}
			else {
				sprintf_s(tipo_inicio, 4, "%s", "AR");
			}
			if (GeneraMasasAzar(MASA_VIA_LACTEA, MASA_VIA_LACTEA) == false) {
				EsperaCaracter();
				return -1;
			}
			if (GeneraCoordenadasAzar(lado_cubo, esfera) == false) {
				EsperaCaracter();
				return -1;
			}
			if (modo_v == 0) {
				if (GeneraVelocidadesHubble(lado_cubo, vel_inicial) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else if (modo_v == 1) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else {
				if (GeneraVelocidadesAzar(true, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (CalculaVelocidadesIntermediasIniciales() == false) {
				EsperaCaracter();
				return -1;
			}
		}
		else if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'R' || ficheromuestra[0] == 'r')) {
			if (modo_v == 0) {
				sprintf_s(tipo_inicio, 4, "%s", "RH");
			}
			else if (modo_v == 1) {
				sprintf_s(tipo_inicio, 4, "%s", "RA");
			}
			else {
				sprintf_s(tipo_inicio, 4, "%s", "RR");
			}
			if (GeneraMasasAzar(MASA_VIA_LACTEA, MASA_VIA_LACTEA) == false) {
				EsperaCaracter();
				return -1;
			}
			if (GeneraCoordenadasReticula(lado_cubo, esfera) == false) {
				EsperaCaracter();
				return -1;
			}
			if (modo_v == 0) {
				if (GeneraVelocidadesHubble(lado_cubo, vel_inicial) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else if (modo_v == 1) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else {
				if (GeneraVelocidadesAzar(true, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (CalculaVelocidadesIntermediasIniciales() == false) {
				EsperaCaracter();
				return -1;
			}
		}
		else {
			int ncd = LeeInicio(ficheromuestra);
			if (ncd == -1) {
				EsperaCaracter();
				return -1;
			}
			if (ncd == 4) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (ncd < 10) {
				if (CalculaVelocidadesIntermediasIniciales() == false) {
					EsperaCaracter();
					return -1;
				}
			}
			sprintf_s(tipo_inicio, 4, "%s", "L_");
		}
	}
	else if (argc == 1) {
		SolicitaCaso();
		SolicitaTipoProceso();
		SolicitaPasoT();
		SolicitaGeometria();
		SolicitaLadoCubo();
		SolicitaVelocidadInicial();
		SolicitaDisColision();
		SolicitaRadioCorte();
		SolicitaFrecuenciaVecinas();
		if (u_vecinas == 1) SolicitaRadioVecinas();
		SolicitaFichero(ficheromuestra);
		SolicitaDuracion();
		SolicitaFrecuenciaMonitor();
		if (t_proceso == 1) {
			SolicitaHilos();
		}
		else {
			n_hilos = 1;
		}
		nseguir = 1;
		particulas_seguir[0] = 0;
		if (!CompletaParametros()) {
			EsperaCaracter();
			return -1;
		}
		if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'A' || ficheromuestra[0] == 'a')) {
			SolicitaGalaxias();
			if (GeneraMasasAzar(MASA_VIA_LACTEA, MASA_VIA_LACTEA) == false) {
				EsperaCaracter();
				return -1;
			}
			if (GeneraCoordenadasAzar(lado_cubo, esfera) == false) {
				return -1;
				EsperaCaracter();
			}
			if (modo_v == 0) {
				if (GeneraVelocidadesHubble(lado_cubo, vel_inicial) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else if (modo_v == 1) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else {
				if (GeneraVelocidadesAzar(true, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (CalculaVelocidadesIntermediasIniciales() == false) {
				EsperaCaracter();
				return -1;
			}
			if (modo_v == 0) {
				sprintf_s(tipo_inicio, 4, "%s", "AH");
			}
			else if (modo_v == 1) {
				sprintf_s(tipo_inicio, 4, "%s", "AA");
			}
			else {
				sprintf_s(tipo_inicio, 4, "%s", "AR");
			}
		}
		else if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'R' || ficheromuestra[0] == 'r')) {
			SolicitaGalaxias();
			if (GeneraMasasAzar(MASA_VIA_LACTEA, MASA_VIA_LACTEA) == false) {
				EsperaCaracter();
				return -1;
			}
			if (GeneraCoordenadasReticula(lado_cubo, esfera) == false) {
				EsperaCaracter();
				return -1;
			}
			if (modo_v == 0) {
				if (GeneraVelocidadesHubble(lado_cubo, vel_inicial) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else if (modo_v == 1) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			else {
				if (GeneraVelocidadesAzar(true, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (CalculaVelocidadesIntermediasIniciales() == false) {
				EsperaCaracter();
				return -1;
			}
			if (modo_v == 0) {
				sprintf_s(tipo_inicio, 4, "%s", "RH");
			}
			else if (modo_v == 1) {
				sprintf_s(tipo_inicio, 4, "%s", "RA");
			}
			else {
				sprintf_s(tipo_inicio, 4, "%s", "RR");
			}
		}
		else {
			int ncd = LeeInicio(ficheromuestra);
			if (ncd == -1) {
				EsperaCaracter();
				return -1;
			}
			if (ncd == 4) {
				if (GeneraVelocidadesAzar(false, vel_inicial, lado_cubo) == false) {
					EsperaCaracter();
					return -1;
				}
			}
			if (ncd < 10) {
				if (CalculaVelocidadesIntermediasIniciales() == false) {
					EsperaCaracter();
					return -1;
				}
			}
			sprintf_s(tipo_inicio, 4, "%s", "L_");
		}
	}
	else {
		printf("Numero de parametros incorrecto, siempre deben ser 15 + particulas a seguir, aunque algunos no se utilicen.\n");
		printf("  numero de caso\n");
		printf("  mono/multi tarea (0 o 1\n");
		printf("  intervalo de tiempo en años\n");
		printf("  geometria (E o C)\n");
		printf("  lado del cubo\n");
		printf("  generar velocidades (velocidad en al/s o H para generarlas con la ley de Hubble, H10 para módulo x10)\n");
		printf("  distancia de colision en al\n");
		printf("  radio de corte en 'al'\n");
		printf("  frecuencia para renovar la tabla de vecinas (0 para no usar tabla)\n");
		printf("  radio de la tabla de vecinas en 'al'\n");
		printf("  datos de la muestra (senda y nombre de fichero, R (en retícula) o A (azar)\n");
		printf("  pasos a simular\n");
		printf("  pasos para monitorizar\n");
		printf("  numero de hilos\n");
		printf("  numero de particulas\n");
		printf("  hasta un máximo de 32, puede no haber ninguna\n");
		printf("	numero de particula a seguir (de 1 hasta numero de particulas)\n");
		printf("	numero de particula a seguir (de 1 hasta numero de particulas)\n");
		printf("	...\n");
		printf("Ejemplos:\n");
		printf("0 0 0.000000031687536 E 1 H 1E-13 2 0 1 F:/Gly/solar.csv 31558150 86400 1 1 1 2 3 4 5 6 7 8 9 10 11 12");
		printf("0 1 1000 E 180000000 H 5E+04 90000000 100 100000000 R 4000000 10000 30 10680");
		EsperaCaracter();
		return -1;
	}
	CANCELADO = false;
	char ficherotmp[256];
	sprintf_s(ficherotmp, 255, filog, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md);
	if ((elog = _fsopen(ficherotmp, "wt", _SH_DENYWR)) == NULL) {
		printf("Error al abrir el fichero log: %s\n", ficherotmp);
		EsperaCaracter();
		return -1;
	}
	if (t_proceso == 0) {
		printf("Inicio proceso mono tarea ...\n");
		ProcesoMono(ficheromuestra, max_iteraciones);
	}
	else {
		printf("Preparando hilos ...\n");
		PreparaMP(ficheromuestra, max_iteraciones);
	}
}

bool CompletaParametros() {
	paso_t *= SEGUNDOS_ANNO;
	medio_paso_t = paso_t / 2.0;
	Gxpaso_t = G * paso_t;
	Gxmedio_paso_t = G * medio_paso_t;
	geometria[0] = esfera ? 'E' : 'C';
	geometria[1] = 0;
	lado_cubo *= AL;
	vel_inicial *= AL;
	dis_critica = dis_colision / 100.0;
	radio_corte *= AL;
	radio_corte2 = radio_corte * radio_corte;
	u_vecinas = f_vecinas == 0 ? 0 : 1;
	radio_vecinas *= AL;
	radio_vecinas2 = radio_vecinas * radio_vecinas;
	frecu_salidas = frecu_monitorizar * 100;
#ifdef UNIDADES_SMI
	sprintf_s(rot_su, 5, "_smi");
#else
	sprintf_s(rot_su, 5, "____");
#endif // UNIDADES_SMI
#ifdef AVX2
	sprintf_s(rot_avx, 6, "_avx2");
#else
	sprintf_s(rot_avx, 6, "_____");
#endif // AVX2
	LiberaMemoria();
	if (!AsignaMemoria(nparticulas)) {
		printf("Error al asignar memoria dinamica");
		return false;
	}
#ifdef MEMORIA_DINAMICA
	sprintf_s(rot_md, 3, "md");
#else
	sprintf_s(rot_md, 3, "__");
#endif // MEMORIA_DINAMICA
	return true;
}

BOOL CtrlHandler(DWORD fdwCtrlType)
{
	// Terminar la ejecución

	CANCELADO = true;

	// Por mucho que pongamos 'Sleep(INFINITE)', el sistema no da más que 10 segundos antes de terminar.
	// Si no se retorna antes de ese tiempo, el sistema cierra de todas formas.

	switch (fdwCtrlType)
	{
	case CTRL_C_EVENT:
		printf("\n\nRecibido CTRL+C. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido CTRL+C. Cerrando ...\n");
		break;
	case CTRL_CLOSE_EVENT:
		printf("\n\nRecibido CTRL+CLOSE. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido CTRL+CLOSE. Cerrando ...\n");
		break;
	case CTRL_BREAK_EVENT:
		printf("\n\nRecibido CTRL+BREAK. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido CTRL+BREAK. Cerrando ...\n");
		break;
	case CTRL_LOGOFF_EVENT:
		printf("\n\nRecibido CTRL+LOGOFF. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido CTRL+LOGOFF. Cerrando ...\n");
		break;
	case CTRL_SHUTDOWN_EVENT:
		printf("\n\nRecibido CTRL+SHUTDOWN. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido CTRL+SHUTDOWN. Cerrando ...\n");
		break;
	default:
		printf("\n\nRecibido Indeterminado. Cerrando, esperar ...\n");
		if (elog != NULL) fprintf(elog, "\nRecibido Indeterminado. Cerrando ...\n");
		break;
	}
	if (elog != NULL) fflush(elog);
	Sleep(INFINITE);
	return FALSE;
}

void EsperaCaracter() {
	printf("\nPulsa retorno de carro\n");
	while (1) {
		if (getchar() == '\n') break;
	}
}

bool BuscaNombreDeProceso(DWORD processID)
{
	// Obtener un Hande al proceso

	HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, processID);

	// Obtener el nombre

	if (NULL != hProcess)
	{
		HMODULE hMod;
		DWORD cbNeeded;

		if (EnumProcessModules(hProcess, &hMod, sizeof(hMod), &cbNeeded)) {
			TCHAR szProcessName[MAX_PATH];
			GetModuleBaseName(hProcess, hMod, szProcessName, sizeof(szProcessName) / sizeof(TCHAR));
			if (_tcscmp(szProcessName, TEXT("Ccababllo.exe")) == 0)
			{
				return true;
			}
		}
	}

	// Liberar el Handle al proceso

	if (hProcess != NULL)CloseHandle(hProcess);
	return false;
}

void SolicitaCaso()
{
	// Pedir el número de caso

	bool correcto;
	char dato[64];
	dato[0] = 0;
	printf("++++ NUMERO DE CASO ++++\n");
	printf(" * = 0\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 5) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				n_caso = 0;;
				break;
			}
			correcto = true;
			for (unsigned int i = 0; i < strlen(dato); i++) {
				if (dato[i] < '0' || dato[i] > '9') {
					printf("Valor incorrecto. Introduce uno nuevo\n");
					correcto = false;
					break;
				}
			}
			if (correcto) {
				n_caso = atoi(dato);
				break;
			}
		}
	}
}

void SolicitaTipoProceso()
{
	// Pedir proceso mono o multitarea

	char dato[64];
	dato[0] = 0;
	printf("++++ TIPO DE PROCESO ++++\n");
	printf(" 0 = Mono\n");
	printf(" 1 = Multi tarea\n");
	printf(" * = Multi tarea.\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) != 1) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (dato[0] == '0') {
				t_proceso = 0;
				break;
			}
			else if (dato[0] == '1') {
				t_proceso = 1;
				break;
			}
			else if (dato[0] == '*') {
				t_proceso = 1;
				break;
			}
			else {
				printf("Valor incorrecto. Introduce uno nuevo\n");
			}
		}
	}
}

void SolicitaPasoT()
{
	// Pedir paso de tiempo en años

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ PASO DE TIEMPO EN ANNOS ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 100 annos\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				paso_t = 1.E+02;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				strcpy_s(datosp, 20, dato);
				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				paso_t = atof(dato);
				if (paso_t < 0.000000001) {
					printf("Valor incorrecto (valor mínimo = 0,000000001). Introduce uno nuevo\n");
				}
				else {
					break;
				}
			}
		}
	}
}

void SolicitaGeometria()
{
	// Pedir geometria esférica o cúbica

	char dato[64];
	dato[0] = 0;
	printf("++++ GEOMETRIA DE LA MUESTRA ++++\n");
	printf(" E = Esterica\n");
	printf(" C = Cubica\n");
	printf(" * = Esferica\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) != 1) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			esfera = dato[0] == '*' || dato[0] == 'E' || dato[0] == 'e';
			break;
		}
	}
}

void SolicitaLadoCubo()
{
	// Pedir el lado del cubo en al

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ LADO DEL CUBO EN AL ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 100.000.000\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				lado_cubo = 1.E+08;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				lado_cubo = atof(dato);
				break;
			}
		}
	}
}

void SolicitaVelocidadInicial() {

	// Pedir la velocidad inicial

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ VELOCIDAD INICIAL al/s ++++\n");
	printf(" H = Hubble, H10 Hubble con el módulo multiplicado por 10\n");
	printf(" Si negativa = Al azar en modulo y direccion radial\n");
	printf(" Se acepta el punto como separador de miles\n");
	printf(" * = 1.0E-10\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				vel_inicial = 1.0E-10;

				// Modulo y dirección al azar

				modo_v = 1;
				break;
			}
			else if (dato[0] == 'h' || dato[0] == 'H') {
				modo_v = 0;
				if (strlen(dato) == 1) {
					vel_inicial = 1.0;
				}
				else {
					dato[0] = '0';
					vel_inicial = atof(ConvPtoComa(dato));
				}
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != '-' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				vel_inicial = atof(dato);
				if (vel_inicial >= 0.0) {

					// Modulo y dirección al azar

					modo_v = 1;
				}
				else {
					// Módulo al azar y dirección radial

					vel_inicial = -vel_inicial;
					modo_v = 2;
				}
				break;
			}
		}
	}
}

void SolicitaDisColision()
{
	// Pedir radio de corte en al

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ DISTANCIA DE COLISION EN al ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 50.000 al\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 16) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {

				// 50 mil al (radio de la via láctea)

				dis_colision = 50000.0;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				dis_colision = atof(dato);
				break;
			}
		}
	}
}

void SolicitaRadioCorte()
{
	// Pedir radio de corte en al

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ RADIO DE CORTE EN AL ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 100.000.000\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				radio_corte = 1.E+08;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				radio_corte = atof(dato);
				break;
			}
		}
	}
}

void SolicitaFrecuenciaVecinas()
{
	// Pedir frecuencia para actualizar la tabla de vecinas

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ FRECUENCIA TABLA DE VECINAS ++++\n");
	printf("  0 = No\n");
	printf(" >0 = Si\n");
	printf("  * = 1000.\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				f_vecinas = 1000;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;
				f_vecinas = atoi(dato);
				if (f_vecinas < 0) {
					printf("Valor incorrecto. Introduce uno nuevo\n");
				}
				else {
					break;
				}
			}
		}
	}
}

void SolicitaRadioVecinas()
{
	// Pedir radio de vecinas en al

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ RADIO VECINAS EN AL ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 100.000.000\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				radio_vecinas = 1.E+08;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != ',' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;

				// Cambiar coma por punto

				n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (dato[i] == ',') {
						dato[i] = '.';
						break;
					}
				}
				radio_vecinas = atof(dato);
				if (radio_vecinas < radio_corte) {
					printf("Valor incorrecto (debe ser mayor o igual que el de corte). Introduce uno nuevo\n");
				}
				else {
					break;
				}
			}
		}
	}
}

char* SolicitaFichero(char* s)
{
	// Pedir fichero de inicio

	char dato[512];
	dato[0] = 0;
	printf(" FICHERO MUESTRA: Masas, Posiciones, Velocidades\n");
	printf(" * = Genera al azar\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (true)
		{
			if (strlen(dato) == 0) {
				printf("Valor incorrecto. Introduce uno nuevo\n");
			}
			else {
				if (strlen(dato) == 1 && dato[0] == '*') {

					// Genera al azar

					dato[0] = 'A';
				}
				break;
			}
		}
	}
	strcpy_s(s, 512, dato);
	return s;
}

void SolicitaGalaxias()
{
	// Pedir número de partículas

	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ NUMERO DE PARTICULAS A GENERAR ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf(" * = 256.\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 12) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				nparticulas = 256;
				break;
			}
			else {
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != '+' && dato[i] != 'E' && dato[i] != 'e') {
						printf("Valor incorrecto. Introduce uno nuevo\n");
						break;
					}
				}

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;
				nparticulas = atoi(dato);
				if (nparticulas > MAX_PARTICULAS) {
					printf("\nDemasiadas partículas (maximo = %d)\n", MAX_PARTICULAS);
				}
				else if (nparticulas < 1) {
					printf("Valor incorrecto. Introduce uno nuevo\n");
				}
				else {
					break;
				}
			}
		}
	}
}

void SolicitaDuracion()
{
	// Pedir duración de la sesión

	bool correcto;
	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ NUMERO MAXIMO DE ITERACIONES PARA ESTA SESION ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf("     Valor en Miles, 1 es mil iteraciones\n");
	printf(" * = 1000 millares  (1 millon)\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 20) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				max_iteraciones = 1000000;
				break;
			}
			correcto = true;
			for (unsigned int i = 0; i < strlen(dato); i++) {
				if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
					printf("Valor incorrecto. Introduce uno nuevo\n");
					correcto = false;
					break;
				}
			}
			if (correcto) {

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;
				max_iteraciones = atol(dato);
				if (max_iteraciones > 0) {
					max_iteraciones *= 1000;
					break;
				}
				else {
					printf("Valor incorrecto. Introduce uno nuevo\n");
				}
			}
		}
	}
}

void SolicitaFrecuenciaMonitor()
{
	// Pedir duración de la sesión

	bool correcto;
	char dato[64];
	char datosp[64];
	dato[0] = 0;
	printf("++++ NUMERO DE ITERACIONES PARA MONITORIZAR ++++\n");
	printf("     Se acepta el punto como separador de miles\n");
	printf("     Valor en Miles, 1 es mil iteraciones\n");
	printf(" * = 100 millares\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 20) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				frecu_monitorizar = 100000;
				break;
			}
			correcto = true;
			for (unsigned int i = 0; i < strlen(dato); i++) {
				if ((dato[i] < '0' || dato[i] > '9') && dato[i] != '.' && dato[i] != '+' && dato[i] != '-' && dato[i] != 'E' && dato[i] != 'e') {
					printf("Valor incorrecto. Introduce uno nuevo\n");
					correcto = false;
					break;
				}
			}
			if (correcto) {

				// Quitar puntos de miles

				strcpy_s(datosp, 20, dato);
				int n = 0;
				for (unsigned int i = 0; i < strlen(dato); i++) {
					if (datosp[i] != '.') dato[n++] = datosp[i];
				}
				dato[n] = 0;
				frecu_monitorizar = atol(dato);
				if (frecu_monitorizar > 0) {
					frecu_monitorizar *= 1000;
					break;
				}
				else {
					printf("Valor incorrecto. Introduce uno nuevo\n");
				}
			}
		}
	}
}

void SolicitaHilos()
{
	// Pedir el número de hilos

	bool correcto;
	char dato[64];
	dato[0] = 0;
	printf("++++ NUMERO DE HILOS ++++\n");
	printf(" * = 20\n");
	printf("Cancelar con Ctrl+C\n");
	while (true) {
		gets_s(dato);
		if (strlen(dato) == 0 || strlen(dato) > 2) {
			printf("Valor incorrecto. Introduce uno nuevo\n");
		}
		else
		{
			if (strlen(dato) == 1 && dato[0] == '*') {
				n_hilos = 20;
				break;
			}
			correcto = true;
			for (unsigned int i = 0; i < strlen(dato); i++) {
				if (dato[i] < '0' || dato[i] > '9') {
					printf("Valor incorrecto. Introduce uno nuevo\n");
					correcto = false;
					break;
				}
			}
			if (correcto) {
				n_hilos = atoi(dato);
				if (n_hilos > 0 && n_hilos <= MAX_HILOS) {
					break;
				}
				else {
					printf("\nDemasiados hilos (maximo = %d)\n", MAX_HILOS);
				}
			}
		}
	}
}

char* ConvPtoComa(char* s) {
	for (int i = 0; i < strlen(s); i++) {
		if (s[i] == '.') s[i] = ',';
	}
	return s;
}

char* SeparaMiles(char* s, unsigned long n) {
	if (n == 0) {
		s[0] = '0';
		s[1] = 0;
		return s;
	}
	unsigned long nori = n;

	unsigned long n2 = 0UL;
	unsigned long scale = 1UL;
	char sp[5];
	s[0] = 0;
	while (n >= 1000UL) {
		n2 = n2 + scale * (n % 1000UL);
		n /= 1000UL;
		scale *= 1000UL;
	}
	sprintf_s(sp, "%lu", n);
	strcat_s(s, 127, sp);
	while (scale != 1) {
		scale /= 1000UL;
		n = n2 / scale;
		n2 = n2 % scale;
		sprintf_s(sp, ".%03lu", n);
		strcat_s(s, 127, sp);
	}
	return s;
}

int DatoColumna(char* linea, int i, char* campo) {

	// Sustituye los '.' por ','

	int uc = (int)strlen(linea) - 1;
	if (i > uc) return -1;
	int j = 0;
	campo[j] = 0;
	while (linea[i] != ';' && i < uc) {
		if (linea[i] == '.') {
			campo[j++] = ',';
		}
		else {
			campo[j++] = linea[i];
		}
		i++;
	}
	campo[j] = 0;
	return i + 1;
}

int LeeInicio(char* fichero) {
	errno_t err;
	FILE* etxt = NULL;
	err = fopen_s(&etxt, fichero, "rt");
	if (err != 0 || etxt == NULL) {
		printf("Error (%d) al abrir el fichero de valores iniciales: %s\n", err, fichero);
		EsperaCaracter();
		return -1;
	}
	char linea[1024];
	char campo[64];

	// Contar el número de líneas, descontando la cabecera

	nparticulas = -1;
	while (feof(etxt) == 0) {
		if (fgets(linea, 1024, etxt) != NULL) nparticulas++;
	}
	if (etxt != NULL)
	{
		if (fclose(etxt))
		{
			printf("No se pudo cerrar el fichero de valores iniciales\n");
			return -1;
		}
	}
	if (nparticulas > MAX_PARTICULAS) {
		printf("\nDemasiadas partículas (maximo = %d)\n", MAX_PARTICULAS);
		EsperaCaracter();
		return -1;
	}
	if (nparticulas < 1) {
		printf("El fichero de valores iniciales esta vacio\n");
		EsperaCaracter();
		return -1;
	}
	err = fopen_s(&etxt, fichero, "rt");
	if (err != 0 || etxt == NULL) {
		printf("Error (%d) al abrir el fichero de valores iniciales: %s\n", err, fichero);
		EsperaCaracter();
		return -1;
	}
	LiberaMemoria();
	if (!AsignaMemoria(nparticulas)) {
		printf("Error al asignar memoria dinamica");
		return -1;
	}

	// Cabecera

	int i;
	int ic = 0;
	int nc = 0;
	fgets(linea, 1024, etxt);
	while (TRUE) {
		if ((ic = DatoColumna(linea, ic, campo)) == -1) {
			break;
		}
		nc++;
	}
	// Quitar la primera columna (nombre)

	nc--;

	if (nc < 4) {
		printf("Numero de columnas insuficiente (%d) debe ser como minimo 4. Error al leer el fichero de valores iniciales: %s\n", nc, fichero);
		EsperaCaracter();
		return -1;
	}

	// Lee los datos

	int n = 0;
	while (feof(etxt) == 0) {
		if (fgets(linea, 1024, etxt) == NULL) continue;
		if (linea == NULL) break;
		i = 0;

		// Columna incrementa el valor de 'i'

		// Ignorar el nombre

		if ((i = DatoColumna(linea, i, campo)) == -1) {
			printf("Error al leer (fila:%d) el campo nombre en el fichero de valores iniciales\n", n);
			EsperaCaracter();
			return -1;
		}

		if ((i = DatoColumna(linea, i, campo)) == -1) {
			printf("Error al leer (fila:%d) el campo masa en el fichero de valores iniciales\n", n);
			EsperaCaracter();
			return -1;
		}
		masas[n] = atof(campo);
		for (int j = 0; j < 3; j++) {
			if ((i = DatoColumna(linea, i, campo)) == -1) {
				printf("Error al leer (fila:%d) el campo coordenada (%d) en el fichero de valores iniciales\n", n, j);
				EsperaCaracter();
				return -1;
			}
			coordenadas[j][n] = atof(campo);
		}
		if (nc > 6) {
			for (int j = 0; j < 3; j++) {
				if ((i = DatoColumna(linea, i, campo)) == -1) {
					printf("Error al leer (fila:%d) el campo velocidad (%d) en el fichero de valores iniciales\n", n, j);
					EsperaCaracter();
					return -1;
				}
				velocidades[j][n] = atof(campo);
			}
		}
		else {
			for (int j = 0; j < 3; j++) {
				velocidades[j][n] = 0.0;
			}
		}
		if (nc > 9) {
			for (int j = 0; j < 3; j++) {
				if ((i = DatoColumna(linea, i, campo)) == -1) {
					printf("Error al leer (fila:%d) el campo velocidad (%d) en el fichero de valores iniciales\n", n, j);
					EsperaCaracter();
					return -1;
				}
				velintermedias[j][n] = atof(campo);
			}
		}
		else {
			for (int j = 0; j < 3; j++) {
				velintermedias[j][n] = 0.0;
			}
		}

		// Las aceleraciones no se leen son función de las masas y coordenadas y por lo tanto se calculan

		n++;
	}
	if (etxt != NULL)
	{
		if (fclose(etxt))
		{
			printf("No se pudo cerrar el fichero de valores iniciales\n");
			return -1;
		}
	}
	return nc;
}

double Covarianza(int n, double x[], double y[])
{
	double xmedia = 0.0;
	double ymedia = 0.0;
	for (int i = 0; i < n; i++)
	{
		xmedia += x[i];
		ymedia += y[i];
	}
	xmedia /= n;
	ymedia /= n;
	double d;
	double xdvstd = 0.0;
	double ydvstd = 0.0;
	double covar = 0.0;
	for (int i = 0; i < n; i++)
	{
		d = (x[i] - xmedia);
		xdvstd += d * d;
		d = (y[i] - ymedia);
		ydvstd += d * d;
		covar += x[i] * y[i];
	}
	xdvstd = sqrt(xdvstd / n);
	ydvstd = sqrt(ydvstd / n);
	covar /= n;
	covar -= xmedia * ymedia;
	covar /= (xdvstd * ydvstd);
	return covar;
}

bool GeneraMasasAzar(double vmin, double vmax) {

	if (vmax - vmin < 1000.0) {
		for (int i = 0; i < nparticulas; i++) {
			masas[i] = vmin;
		}
		return true;
	}
	mt19937 mt(1235);
	uniform_real_distribution<double> dist(vmin, vmax);
	double v;
	for (int i = 0; i < nparticulas; i++) {
		v = dist(mt);
		masas[i] = v;
	}
	return true;
}

bool GeneraCoordenadasAzar(double lcubo, bool esfera) {

	// Centradas en 0, 0, 0

	double v;
	double d;
	double d2;
	int ind;
	double cc[3];
	double radio = lcubo / 2.0;
	mt19937 mt(17555);
	uniform_real_distribution<double> distribucion(0.0, 1.0);
	const int NUM_HISTO_D = NUM_HISTO;
	double dis[NUM_HISTO_D];
	double vol[NUM_HISTO_D];
	int histo_d[NUM_HISTO_D];
	int p_exactas[NUM_HISTO_D];
	double delta_r = radio / (double)NUM_HISTO_D;
	int max_intentos = 0;
	if (esfera) {
		double pi43 = 4.0 / 3.0 * M_PI;
		double volumen = pi43 * radio * radio * radio;
		double densidad = (double)nparticulas / volumen;
		double difabs;
		double r1 = 0.0;
		double r2;
		for (int i = 0; i < NUM_HISTO_D; i++) {
			r2 = r1 + delta_r;
			dis[i] = (r1 + r2) / 2.0;
			vol[i] = pi43 * (r2 * r2 * r2 - r1 * r1 * r1);
			p_exactas[i] = (int)round(densidad * vol[i]);
			r1 += delta_r;
			histo_d[i] = 0;
		}
		int np;
		bool seguir = true;
		while (seguir == true && max_intentos < 1000) {
			max_intentos++;
			np = 0;
			while (np < nparticulas) {
				d2 = 0.0;
				for (int k = 0; k < 3; k++) {
					v = distribucion(mt) * lcubo - radio;
					cc[k] = v;
					d2 += v * v;
				}
				d = sqrt(d2);
				if (d <= radio) {
					for (int k = 0; k < 3; k++) {
						coordenadas[k][np] = cc[k];
					}
					np++;
				}
			}
			for (int i = 0; i < NUM_HISTO_D; i++) {
				histo_d[i] = 0;
			}
			for (int i = 0; i < nparticulas; i++) {
				d2 = 0.0;
				for (int k = 0; k < 3; k++) {
					d2 += coordenadas[k][i] * coordenadas[k][i];
				}
				d = sqrt(d2);
				ind = (int)(d / delta_r);
				histo_d[ind]++;
			}
			seguir = false;
			for (int i = 0; i < NUM_HISTO_D; i++) {
				difabs = abs(histo_d[i] - p_exactas[i]);
				if (difabs / (double)p_exactas[i] > 0.1) {
					seguir = true;
					break;
				}
			}
		}
	}
	else {
		for (int i = 0; i < nparticulas; i++) {
			for (int k = 0; k < 3; k++) {
				v = distribucion(mt) * lcubo - radio;
				coordenadas[k][i] = v;
			}
		}
	}

	// Para comprobar la uniformidad de la muestra

	double minimo[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
	double maximo[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
	double media[3] = { 0.0,0.0,0.0 };
	for (int i = 0; i < nparticulas; i++) {
		for (int k = 0; k < 3; k++) {
			if (coordenadas[k][i] < minimo[k]) minimo[k] = coordenadas[k][i];
			if (coordenadas[k][i] > maximo[k]) maximo[k] = coordenadas[k][i];
			media[k] += coordenadas[k][i];
		}
	}
	for (int k = 0; k < 3; k++) {
		minimo[k] /= radio;
		maximo[k] /= radio;
	}
	for (int k = 0; k < 3; k++) {
		media[k] /= nparticulas;
	}
	double covar[3];
	covar[0] = Covarianza(nparticulas, coordenadas[0], coordenadas[1]);
	covar[1] = Covarianza(nparticulas, coordenadas[0], coordenadas[2]);
	covar[2] = Covarianza(nparticulas, coordenadas[1], coordenadas[2]);
	double sum[3] = { 0.0,0.0,0.0 };
	for (int i = 0; i < nparticulas; i++) {
		for (int k = 0; k < 3; k++) {
			v = coordenadas[k][i] - media[k];
			sum[k] += v * v;
		}
	}
	double de[3];
	for (int k = 0; k < 3; k++) {
		de[k] = sqrt(sum[k] / ((double)nparticulas - 1));
	}
	if (NUM_HISTO > 0) {
		int histo[NUM_HISTO][3];
		double delta = lcubo / NUM_HISTO;
		for (int i = 0; i < NUM_HISTO; i++) {
			for (int k = 0; k < 3; k++) {
				histo[i][k] = 0;
			}
		}
		for (int i = 0; i < nparticulas; i++) {
			for (int k = 0; k < 3; k++) {
				ind = (int)((coordenadas[k][i] + radio) / delta);
				histo[ind][k]++;
			}
		}
		errno_t err;
		FILE* ec;
		char ficherotmp[256];
		sprintf_s(ficherotmp, 255, fihtc, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md);
		err = fopen_s(&ec, ficherotmp, "wt");
		if (err != 0 || ec == NULL) {
			printf("Error (%d) al abrir el fichero para el histograma de coordenadas: %s\n", err, ficherotmp);
			EsperaCaracter();
			return false;
		}
		fprintf(ec, "x;y;z\n");
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", minimo[0], minimo[1], minimo[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", maximo[0], maximo[1], maximo[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", media[0] / lcubo, media[1] / lcubo, media[2] / lcubo);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", de[0] / lcubo, de[1] / lcubo, de[2] / lcubo);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", covar[0], covar[1], covar[2]);
		for (int i = 0; i < NUM_HISTO; i++) {
			fprintf(ec, "%5d;%5d;%5d\n", histo[i][0], histo[i][1], histo[i][2]);
		}
		if (esfera) {
			for (int i = 0; i < NUM_HISTO_D; i++) {
				histo_d[i] = 0;
			}
			for (int i = 0; i < nparticulas; i++) {
				d2 = 0.0;
				for (int k = 0; k < 3; k++) {
					d2 += coordenadas[k][i] * coordenadas[k][i];
				}
				d = sqrt(d2);
				ind = (int)(d / delta_r);
				histo_d[ind]++;
			}
			fprintf(ec, "\ndistancia;volumen;teoricas;particulas;densidad\n");
			for (int i = 0; i < NUM_HISTO_D; i++) {
				fprintf(ec, "%14.7e;%14.7e;%5d;%5d;%14.7e\n", dis[i], vol[i], p_exactas[i], histo_d[i], (double)histo_d[i] / vol[i]);
			}
		}
		if (ec != NULL) fclose(ec);
	}
	return true;
}

bool GeneraCoordenadasReticula(double lcubo, bool esfera) {

	// Centradas en 0, 0, 0

	double radio = lcubo / 2.0;
	double x;
	double y;
	double z;
	double d;
	int n;
	int np = 0;
	double delta;
	if (esfera) {
		double rce = (lcubo * lcubo * lcubo) / (4.0 / 3.0 * M_PI * (radio * radio * radio));
		n = (int)pow(nparticulas * rce, 1.0 / 3.0) - 1;
		while (np < nparticulas) {
			n++;
			np = 0;
			delta = lcubo / (double)n;
			for (int i = 0; i < n; i++) {
				x = (i + 0.5) * delta - radio;
				for (int j = 0; j < n; j++) {
					y = (j + 0.5) * delta - radio;
					for (int k = 0; k < n; k++) {
						z = (k + 0.5) * delta - radio;
						d = sqrt(x * x + y * y + z * z);
						if (d <= radio) {
							np++;
							if (np > nparticulas) break;
						}
					}
					if (np > nparticulas) break;
				}
				if (np > nparticulas) break;
			}
		}
		if (np > nparticulas) n--;
	}
	else {
		n = (int)pow(nparticulas, 1.0 / 3.0);
	}
	if (n < 2) {
		printf("\nNumero de particulas insuficiente.\n");
		return false;
	}
	np = 0;
	delta = lcubo / (double)n;
	for (int i = 0; i < n; i++) {
		x = (i + 0.5) * delta - radio;
		for (int j = 0; j < n; j++) {
			y = (j + 0.5) * delta - radio;
			for (int k = 0; k < n; k++) {
				z = (k + 0.5) * delta - radio;
				if (esfera) {
					d = sqrt(x * x + y * y + z * z);
					if (d <= radio) {
						coordenadas[0][np] = x;
						coordenadas[1][np] = y;
						coordenadas[2][np] = z;
						np++;
						if (np == nparticulas) break;
					}
				}
				else {
					coordenadas[0][np] = x;
					coordenadas[1][np] = y;
					coordenadas[2][np] = z;
					np++;
				}
			}
			if (np == nparticulas) break;
		}
		if (np == nparticulas) break;
	}
	nficticias = nparticulas - np;
	if (nficticias > 0) {

		// Poner el resto en el origen sin masa

		for (int i = np; i < nparticulas; i++) {
			masas[i] = 0.0;
			for (int k = 0; k < 3; k++) {
				coordenadas[k][i] = 0.0;
			}
		}
	}
	return true;
}

bool GeneraVelocidadesAzar(bool radial, double vmax, double lcubo) {

	// radial=true :En la dirección del vector de posición 'r' respecto al origen

	mt19937 mt(3512);
	uniform_real_distribution<double> dist(0.0, 1.0);
	double ur2;
	double ur;
	double mod_v;
	double vuv[3];
	double uv2;
	double uv;
	double ve;
	double minimo[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
	double maximo[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
	double media[3] = { 0.0,0.0,0.0 };
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) {
			for (int k = 0; k < 3; k++) {
				velocidades[k][i] = 0.0;
			}
		}
		else {
			if (radial) {
				ur2 = 0.0;
				for (int k = 0; k < 3; k++) {
					ur2 += coordenadas[k][i] * coordenadas[k][i];
				}
				ur = sqrt(ur2);
				if (ur < 1.0E-20) {
					for (int k = 0; k < 3; k++) {
						mod_v = 0.0;
						if (mod_v < minimo[k]) minimo[k] = mod_v;
						if (mod_v > maximo[k]) maximo[k] = mod_v;
						velocidades[k][i] = mod_v;
					}
				}
				else {
					for (int k = 0; k < 3; k++) {
						mod_v = dist(mt) * vmax;
						ve = mod_v * coordenadas[k][i] / ur;
						velocidades[k][i] = ve;
						if (mod_v < minimo[k]) minimo[k] = mod_v;
						if (mod_v > maximo[k]) maximo[k] = mod_v;
						media[k] += ve;
					}
				}
			}
			else {
				for (int k = 0; k < 3; k++) {

					// Dirección al azar

					uv2 = 0.0;
					for (int k = 0; k < 3; k++) {
						vuv[k] = dist(mt);
						uv2 += vuv[k] * vuv[k];
					}
					uv = sqrt(uv2);

					// Módulo

					mod_v = dist(mt) * vmax;
					ve = mod_v * vuv[k] / uv;
					velocidades[k][i] = ve;
					if (mod_v < minimo[k]) minimo[k] = mod_v;
					if (mod_v > maximo[k]) maximo[k] = mod_v;
					media[k] += ve;
				}
			}
		}
	}

	// Para comprobar la uniformidad de la muestra

	for (int k = 0; k < 3; k++) {
		minimo[k] /= vmax;
		maximo[k] /= vmax;
	}
	for (int k = 0; k < 3; k++) {
		media[k] /= nparticulas;
	}
	double covar[3];
	covar[0] = Covarianza(nparticulas, velocidades[0], velocidades[1]);
	covar[1] = Covarianza(nparticulas, velocidades[0], velocidades[2]);
	covar[2] = Covarianza(nparticulas, velocidades[1], velocidades[2]);
	double sum[3] = { 0.0,0.0,0.0 };
	for (int i = 0; i < nparticulas; i++) {
		for (int k = 0; k < 3; k++) {
			mod_v = velocidades[k][i] - media[k];
			sum[k] += mod_v * mod_v;
		}
	}
	double de[3];
	for (int k = 0; k < 3; k++) {
		de[k] = sqrt(sum[k] / ((double)nparticulas - 1));
	}
	if (NUM_HISTO > 0) {
		int histo[NUM_HISTO * 2][3];
		for (int i = 0; i < NUM_HISTO * 2; i++) {
			for (int k = 0; k < 3; k++) {
				histo[i][k] = 0;
			}
		}
		double delta = vmax / NUM_HISTO;
		int ind;
		for (int i = 0; i < nparticulas; i++) {
			for (int k = 0; k < 3; k++) {
				ind = (int)(velocidades[k][i] / delta) + NUM_HISTO;
				histo[ind][k]++;
			}
		}
		errno_t err;
		FILE* ec;
		char ficherotmp[256];
		sprintf_s(ficherotmp, 255, fihtv, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md);
		err = fopen_s(&ec, ficherotmp, "wt");
		if (err != 0 || ec == NULL) {
			printf("Error (%d) al abrir el fichero para el histograma de velocidades: %s\n", err, ficherotmp);
			EsperaCaracter();
			return false;
		}
		fprintf(ec, "x;y;z\n");
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", minimo[0], minimo[1], minimo[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", maximo[0], maximo[1], maximo[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", media[0], media[1], media[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", de[0], de[1], de[2]);
		fprintf(ec, "%14.7e;%14.7e;%14.7e\n", covar[0], covar[1], covar[2]);
		for (int i = 0; i < NUM_HISTO * 2; i++) {
			fprintf(ec, "%5d;%5d;%5d\n", histo[i][0], histo[i][1], histo[i][2]);
		}
		if (ec != NULL) fclose(ec);
	}
	return true;
}

bool GeneraVelocidadesHubble(double lcubo, double vel_inicial) {

	// En la dirección del vector de posición 'r' respecto del origen

	double d2;
	double d;
	double v;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) {
			for (int k = 0; k < 3; k++) {
				velocidades[k][i] = 0.0;
			}
		}
		else {

			// Distancia al centro y dirección de la velocidad

			d2 = 0.0;
			for (int k = 0; k < 3; k++) {
				d2 += coordenadas[k][i] * coordenadas[k][i];
			}
			d = sqrt(d2);

			// Módulo y componentes de la velocidad

			if (d < 1.0E-20) {
				for (int k = 0; k < 3; k++) {
					velocidades[k][i] = 0.0;
				}
			}
			else {

				// Hubble: v = H0 * D 

				v = H0 * d;
				for (int k = 0; k < 3; k++) {
					velocidades[k][i] = vel_inicial * v * coordenadas[k][i] / d;
				}
			}
		}
	}
	return true;
}

bool CalculaVelocidadesIntermediasIniciales() {

	// Calcular las velocidades intermedias de partida

	Aceleraciones();
	if (nparticulas == 1) {
		printf("El universo está colapsado\n");
		return false;
	}
	for (int n = 0; n < nparticulas; n++) {
		for (int k = 0; k < 3; k++) {
			velintermedias[k][n] = velocidades[k][n] - aceleraciones[k][n] * Gxmedio_paso_t;
		}
	}
	return true;
}

bool EscribeDatos(int n) {
	errno_t err;
	FILE* ec;
	char ficherotmp[256];
	sprintf_s(ficherotmp, 255, fisal, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md, n);
	err = fopen_s(&ec, ficherotmp, "wt");
	if (err != 0 || ec == NULL) {
		printf("Error (%d) al abrir el fichero de salida: %s\n", err, ficherotmp);
		EsperaCaracter();
		return false;
	}
	fprintf(ec, "id;m;x;y;z;vx;vy;vz;vix;viy;viz;ax;ay;az\n");
	for (int i = 0; i < nparticulas; i++) {
		ImprimeDatosParticula(ec, i);
	}
	if (ec != NULL) fclose(ec);
	return true;
}

double MediaVecinas() {
	if (t_vecinas_enuso == -1) return 0.0;
	int n = 0;
	double media = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		media += (double)n_vecinas[i][t_vecinas_enuso];
		n++;
	}
	if (n == 0) {
		return 0.0;
	}
	else {
		return media / (double)n;
	}
}

void CentroMasas(double* cm) {

	// Centro de masas

	double sm = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		cm[0] += masas[i] * coordenadas[0][i];
		cm[1] += masas[i] * coordenadas[1][i];
		cm[2] += masas[i] * coordenadas[2][i];
		sm += masas[i];
	}
	cm[0] /= sm;
	cm[1] /= sm;
	cm[2] /= sm;
}

double DistanciaMediaAcm(double* cm) {

	// Distancia mdia ponderada respecto al centro de masas

	double x;
	double y;
	double z;
	double dm = 0.0;
	double sm = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		x = coordenadas[0][i] - cm[0];
		y = coordenadas[1][i] - cm[1];
		z = coordenadas[2][i] - cm[2];
		dm += masas[i] * sqrt(x * x + y * y + z * z);
		sm += masas[i];
	}
	return dm / sm;
}

double DesviacionEstandarDisAcm(double* cm, double dm) {

	// Desviación estándar

	double x;
	double y;
	double z;
	double dd;
	double de = 0.0;
	double sm = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		x = coordenadas[0][i] - cm[0];
		y = coordenadas[1][i] - cm[1];
		z = coordenadas[2][i] - cm[2];
		dd = sqrt(x * x + y * y + z * z) - dm;
		de += masas[i] * dd * dd;
		sm += masas[i];
	}
	return sqrt(de / sm);
}

double Ecinetica() {
	double v2;
	double ec = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		v2 = velocidades[0][i] * velocidades[0][i] + velocidades[1][i] * velocidades[1][i] + velocidades[2][i] * velocidades[2][i];
		ec += masas[i] * v2 / 2.0;
	}
	return ec;
}

double EcineticaParticula(int i) {
	if (masas[i] == 0.0) return 0.0;
	double v2 = velocidades[0][i] * velocidades[0][i] + velocidades[1][i] * velocidades[1][i] + velocidades[2][i] * velocidades[2][i];
	return masas[i] * v2 / 2.0;
}

double Epotencial(double* epr, double* ndentro) {
	double x;
	double y;
	double z;
	double d2;
	double d;
	double val;
	double ep = 0.0;
	int n_d = 0;
	int n_p = 0;
	*epr = 0.0;
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		for (int j = i + 1; j < nparticulas; j++) {
			if (masas[j] == 0.0) break;
			x = coordenadas[0][i] - coordenadas[0][j];
			y = coordenadas[1][i] - coordenadas[1][j];
			z = coordenadas[2][i] - coordenadas[2][j];
			d2 = x * x + y * y + z * z;
			d = sqrt(d2);
			if (d < dis_critica) {
				val = -G * masas[i] * masas[j] / dis_critica;
			}
			else {
				val = -G * masas[i] * masas[j] / d;
			}
			if (d > radio_corte) {
				(*epr) += val;
			}
			else {
				ep += val;
				n_d++;
			}
		}
		n_p++;
	}
	if (n_p == 0) {
		(*ndentro) = 0.0;
	}
	else {
		(*ndentro) = (double)n_d * 2.0 / (double)n_p;
	}
	return ep;
}

double EpotencialParticula(int i, double* epr) {
	*epr = 0.0;
	if (masas[i] == 0.0) return 0.0;
	double x;
	double y;
	double z;
	double d2;
	double d;
	double val;
	double ep = 0.0;
	for (int j = i + 1; j < nparticulas; j++) {
		if (masas[j] == 0.0) break;
		x = coordenadas[0][i] - coordenadas[0][j];
		y = coordenadas[1][i] - coordenadas[1][j];
		z = coordenadas[2][i] - coordenadas[2][j];
		d2 = x * x + y * y + z * z;
		d = sqrt(d2);
		if (d < dis_critica) {
			val = -G * masas[i] * masas[j] / dis_critica;
		}
		else {
			val = -G * masas[i] * masas[j] / d;
		}
		if (d > radio_corte) {
			(*epr) += val;
		}
		else {
			ep += val;
		}
	}
	return ep;
}

bool Arranque(char* ficheromuestra, unsigned long hasta) {
	time(&reloji);
	ncolisiones = 0;
	ncriticas = 0;

	// Escribir los datos de partida

	nsalida = 0;
	EscribeDatos(nsalida++);
#ifdef AVX2
	fprintf(elog, "Simulación Movimiento de partículas Gravitatorias. v %d.%02d AVX 256\n", version, subversion);
#else
	fprintf(elog, "Simulación Movimiento de partículas Gravitatorias. v %d.%02d\n", version, subversion);
#endif // AVX2
	fprintf(elog, "Partículas          %14d\n", nparticulas);
	if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'A' || ficheromuestra[0] == 'a')) {
		fprintf(elog, "Inicio al azar\n");
	}
	else if (strlen(ficheromuestra) == 1 && (ficheromuestra[0] == 'R' || ficheromuestra[0] == 'r')) {
		fprintf(elog, "Inicio desde retícula. %d ficticias\n", nficticias);
	}
	else {
		fprintf(elog, "Inicio desde %s\n", ficheromuestra);
	}
	if (t_proceso == 0) {
		fprintf(elog, "Mono taréa\n");
	}
	else {
		fprintf(elog, "Número de hilos     %14d\n", n_hilos);
	}
	fprintf(elog, "Paso de tiempo      %14.7e años\n", paso_t / SEGUNDOS_ANNO);
	if (esfera) {
		fprintf(elog, "Muestra esférica\n");
		fprintf(elog, "Radio de la esfera  %14.7e   al\n", lado_cubo / (2.0 * AL));
	}
	else {
		fprintf(elog, "Muestra cúbica\n");
		fprintf(elog, "Lado del cubo       %14.7e   al\n", lado_cubo / AL);
	}
	if (modo_v == 0) {
		fprintf(elog, "Velocidades Hubble  %14.7e factor\n", vel_inicial);
	}
	else if (modo_v == 1) {
		fprintf(elog, "Velocidad máxima    %14.7e al/s al azar en módulo y dirección\n", vel_inicial);
	}
	else {
		fprintf(elog, "Velocidad máxima    %14.7e al/s al azar en módulo y dirección radial\n", vel_inicial);
	}
	fprintf(elog, "Distancia colisión  %14.7e   al\n", dis_colision / AL);
	fprintf(elog, "Radio de corte      %14.7e   al\n", radio_corte / AL);
	if (u_vecinas == 0) {
		fprintf(elog, "Tabla de vecinas    %14s\n", "NO");
	}
	else {
#ifdef HILO_VECINAS
		fprintf(elog, "Radio vecinas       %14.7e   al\n", radio_vecinas / AL);
		fprintf(elog, "Frecuencia vecinas  %14s Hilo\n", SeparaMiles(fto1, f_vecinas));
		ActualizaTablaVecinas();
#else
		fprintf(elog, "Frecuencia vecinas  %14s\n", SeparaMiles(fto1, f_vecinas));
#endif //  HILO_VECINAS
	}
	fprintf(elog, "Monitorizar cada    %14s pasos\n", SeparaMiles(fto1, frecu_monitorizar));
	fprintf(elog, "Tiempo a simular    %14.8f millones de años (%s pasos)\n\n", hasta * paso_t / SEGUNDOS_ANNO / 1000000.0, SeparaMiles(fto1, hasta));

	char ficherotmp[256];

	// Fichero de colisiones

	sprintf_s(ficherotmp, 255, ficol, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md);
	if ((ecol = _fsopen(ficherotmp, "wt", _SH_DENYWR)) == NULL) {
		printf("Error al abrir el fichero de colisiones: %s\n", ficherotmp);
		EsperaCaracter();
		return false;
	}

	// Seguimiento de 'nseguir' partículas

	int j;
	for (int i = 0; i < nseguir; i++) {
		j = particulas_seguir[i];
		coor_ini[0][i] = coordenadas[0][j];
		coor_ini[1][i] = coordenadas[1][j];
		coor_ini[2][i] = coordenadas[2][j];
		recorrido[i] = 0.0;
		sprintf_s(ficherotmp, 255, fiseg, rot_su, rot_avx, t_proceso, n_hilos, tipo_inicio, geometria, nparticulas, n_caso, rot_md, i + 1);
		if ((eseg[i] = _fsopen(ficherotmp, "wt", _SH_DENYWR)) == NULL) {
			printf("Error al abrir el fichero de seguimiento: %s\n", ficherotmp);
			EsperaCaracter();
			return false;
		}
		fprintf(eseg[i], "t;m;x;y;z;vx;vy;vz;vix;viy;viz;ax;ay;az;v;ec;ep;epr;rec;des;n;par\n");
	}
	ImprimeSeguimiento(0);

	double v_pp = sqrt(velocidades[0][0] * velocidades[0][0] + velocidades[1][0] * velocidades[1][0] + velocidades[2][0] * velocidades[2][0]);
	double ec_pp = EcineticaParticula(0);
	double epr_pp;
	double ep_pp = EpotencialParticula(0, &epr_pp);
	double et_pp = ec_pp + ep_pp;
	fprintf(elog, "                ener.cinetica ener.potencial  ENERGIA TOTAL ener.pot.resid           x al           y al           z al        vx al/s        vy al/s        vz al/s         v al/s\n");
	fprintf(elog, "-------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	fprintf(elog, "Partícula uno  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", ec_pp, ep_pp, et_pp, epr_pp, coordenadas[0][0] / AL, coordenadas[1][0] / AL, coordenadas[2][0] / AL, velocidades[0][0] / AL, velocidades[1][0] / AL, velocidades[2][0] / AL, v_pp / AL);
	fprintf(elog, "-------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n\n");

	printf(" millones anos  ener.cinetica ener.potencial  ENERGIA TOTAL In.corte ener.pot.resid   recorrido al desplazami. al      X c.m. al      Y c.m. al      Z c.m. al dist. media al desv.estand al  tpo calculo s\n");
	printf("-------------- -------------- -------------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	fprintf(elog, " millones años  ener.cinética ener.potencial  ENERGÍA TOTAL In.corte ener.pot.resid In.veci.   recorrido al desplazami. al      X c.m. al      Y c.m. al      Z c.m. al dist. media al desv.estand al  tpo cálculo s\n");
	fprintf(elog, "-------------- -------------- -------------- -------------- -------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	double nvemed = MediaVecinas();
	j = particulas_seguir[0];
	double dx = coordenadas[0][j] - coor_ini[0][0];
	double dy = coordenadas[1][j] - coor_ini[1][0];
	double dz = coordenadas[2][j] - coor_ini[2][0];
	double desplazamiento = sqrt(dx * dx + dy * dy + dz * dz);
	double cm[3] = { 0.0,0.0,0.0 };
	double ec = Ecinetica();
	double epr;
	double ndentro;
	double ep = Epotencial(&epr, &ndentro);
	double et = ec + ep;
	CentroMasas(cm);
	double dm = DistanciaMediaAcm(cm);
	double de = DesviacionEstandarDisAcm(cm, dm);
	time(&reloja);
	segundos = difftime(reloja, reloji);
	printf("%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	fprintf(elog, "%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, nvemed, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	fflush(elog);
	return true;
}

void Seguimiento(unsigned long iteracion) {
	double nvemed = MediaVecinas();
	int j = particulas_seguir[0];
	double dx = coordenadas[0][j] - coor_ini[0][0];
	double dy = coordenadas[1][j] - coor_ini[1][0];
	double dz = coordenadas[2][j] - coor_ini[2][0];
	double desplazamiento = sqrt(dx * dx + dy * dy + dz * dz);
	double ec = Ecinetica();
	double epr;
	double ndentro;
	double ep = Epotencial(&epr, &ndentro);
	double et = ec + ep;
	double cm[3] = { 0.0,0.0,0.0 };
	CentroMasas(cm);
	double dm = DistanciaMediaAcm(cm);
	double de = DesviacionEstandarDisAcm(cm, dm);
	time(&reloja);
	segundos = difftime(reloja, reloji);
	printf("%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	fprintf(elog, "%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, nvemed, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	fflush(elog);
	ImprimeSeguimiento(iteracion);
}

void Terminacion(unsigned long iteracion, unsigned long hasta, const char* rotulo) {
	double nvemed = MediaVecinas();
	int j = particulas_seguir[0];
	double dx = coordenadas[0][j] - coor_ini[0][0];
	double dy = coordenadas[1][j] - coor_ini[1][0];
	double dz = coordenadas[2][j] - coor_ini[2][0];
	double desplazamiento = sqrt(dx * dx + dy * dy + dz * dz);
	double ec = Ecinetica();
	double epr;
	double ndentro;
	double ep = Epotencial(&epr, &ndentro);
	double et = ec + ep;
	double cm[3] = { 0.0,0.0,0.0 };
	CentroMasas(cm);
	double dm = DistanciaMediaAcm(cm);
	double de = DesviacionEstandarDisAcm(cm, dm);
	time(&reloja);
	segundos = difftime(reloja, reloji);
	printf("-------------- -------------- -------------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	printf("%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	printf("-------------- -------------- -------------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	fprintf(elog, "-------------- -------------- -------------- -------------- -------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	fprintf(elog, "%14.8f %14.7e %14.7e %14.7e %8.2f %14.7e %8.2f %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", iteracion * paso_t / SEGUNDOS_ANNO / 1000000.0, ec, ep, et, ndentro, epr, nvemed, recorrido[0] / AL, desplazamiento / AL, cm[0] / AL, cm[1] / AL, cm[2] / AL, dm / AL, de / AL, segundos);
	fprintf(elog, "-------------- -------------- -------------- -------------- -------- -------------- -------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	double v_pp = sqrt(velocidades[0][0] * velocidades[0][0] + velocidades[1][0] * velocidades[1][0] + velocidades[2][0] * velocidades[2][0]);
	double ec_pp = EcineticaParticula(0);
	double epr_pp;
	double ep_pp = EpotencialParticula(0, &epr_pp);
	double et_pp = ec_pp + ep_pp;
	fprintf(elog, "                ener.cinetica ener.potencial  ENERGIA TOTAL ener.pot.resid           x al           y al           z al        vx al/s        vy al/s        vz al/s         v al/s\n");
	fprintf(elog, "-------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");
	fprintf(elog, "Partícula uno  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", ec_pp, ep_pp, et_pp, epr_pp, coordenadas[0][0] / AL, coordenadas[1][0] / AL, coordenadas[2][0] / AL, velocidades[0][0] / AL, velocidades[1][0] / AL, velocidades[2][0] / AL, v_pp / AL);
	fprintf(elog, "-------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------\n");

	// Escribir los datos finales

	EscribeDatos(nsalida);
	time(&reloja);
	segundos = difftime(reloja, reloji);
	printf("\n%d colisiones, %d criticas \n", ncolisiones, ncriticas);
	printf("Numero maximo de vecinas %d\n", max_absoluto_vecinas);
	printf("\n%s. %s pasos.  Tiempo: %s s\n", rotulo, SeparaMiles(fto1, iteracion), SeparaMiles(fto2, (unsigned long)segundos));
	fprintf(elog, "\n%d colisiones, %d críticas\n", ncolisiones, ncriticas);
	fprintf(elog, "Número máximo de vecinas %d\n", max_absoluto_vecinas);
	fprintf(elog, "\n%s. %s pasos.  Tiempo: %s s\n", rotulo, SeparaMiles(fto1, iteracion), SeparaMiles(fto2, (unsigned long)segundos));
	Beep(800, 500);
	fclose(ecol);
	fclose(elog);
	for (int i = 0; i < nseguir; i++) {
		fclose(eseg[i]);
	}
}

void MuestraDatosParticula(FILE* ef, int i) {
	fprintf(ef, "%14.7e;", masas[i]);
	fprintf(ef, "%14.7e;%14.7e;%14.7e;", coordenadas[0][i], coordenadas[1][i], coordenadas[2][i]);
	fprintf(ef, "%14.7e;%14.7e;%14.7e;", velocidades[0][i], velocidades[1][i], velocidades[2][i]);
	fprintf(ef, "%14.7e;%14.7e;%14.7e;", velintermedias[0][i], velintermedias[1][i], velintermedias[2][i]);
	fprintf(ef, "%14.7e;%14.7e;%14.7e\n", aceleraciones[0][i], aceleraciones[1][i], aceleraciones[2][i]);
}

void ImprimeDatosParticula(FILE* ef, int i) {
	fprintf(ef, "%16d;", i);
	fprintf(ef, "%16.9e;", masas[i]);
	fprintf(ef, "%16.9e;%16.9e;%16.9e;", coordenadas[0][i], coordenadas[1][i], coordenadas[2][i]);
	fprintf(ef, "%16.9e;%16.9e;%16.9e;", velocidades[0][i], velocidades[1][i], velocidades[2][i]);
	fprintf(ef, "%16.9e;%16.9e;%16.9e;", velintermedias[0][i], velintermedias[1][i], velintermedias[2][i]);
	fprintf(ef, "%16.9e;%16.9e;%16.9e\n", aceleraciones[0][i], aceleraciones[1][i], aceleraciones[2][i]);
}

void ImprimeSeguimiento(unsigned long iteracion) {
	double v_pp;
	double ec_pp;
	double epr_pp;
	double ep_pp;
	double desplazamiento;
	for (int i = 0; i < nseguir; i++) {
		int j = particulas_seguir[i];
		double dx = coordenadas[0][j] - coor_ini[0][i];
		double dy = coordenadas[1][j] - coor_ini[1][i];
		double dz = coordenadas[2][j] - coor_ini[2][i];
		desplazamiento = sqrt(dx * dx + dy * dy + dz * dz);
		v_pp = sqrt(velocidades[0][j] * velocidades[0][j] + velocidades[1][j] * velocidades[1][j] + velocidades[2][j] * velocidades[2][j]);
		ec_pp = EcineticaParticula(j);
		ep_pp = EpotencialParticula(j, &epr_pp);
		fprintf(eseg[i], "%16.9e;%16.9e;", iteracion * paso_t, masas[i]);
		fprintf(eseg[i], "%16.9e;%16.9e;%16.9e;", coordenadas[0][j], coordenadas[1][j], coordenadas[2][j]);
		fprintf(eseg[i], "%16.9e;%16.9e;%16.9e;", velocidades[0][j], velocidades[1][j], velocidades[2][j]);
		fprintf(eseg[i], "%16.9e;%16.9e;%16.9e;", velintermedias[0][j], velintermedias[1][j], velintermedias[2][j]);
		fprintf(eseg[i], "%16.9e;%16.9e;%16.9e;", aceleraciones[0][j], aceleraciones[1][j], aceleraciones[2][j]);
		fprintf(eseg[i], "%16.9e;%16.9e;%16.9e;%16.9e;", v_pp, ec_pp, ep_pp, epr_pp);
		fprintf(eseg[i], "%16.9e;%16.9e;%3d;%6d\n", recorrido[i], desplazamiento, i + 1, j + 1);
		fflush(eseg[i]);
	}
}

void ProcesoMono(char* ficheromuestra, unsigned long hasta) {
	if (!Arranque(ficheromuestra, hasta)) return;
#ifdef AVX2
	npt_mt = nparticulas;
	npm4_mt = (npt_mt / bloque_mt) * bloque_mt;
	pasot_mt = _mm256_set1_pd(paso_t);
	gxpasot_mt = _mm256_set1_pd(Gxpaso_t);
	gxmediopasot_mt = _mm256_set1_pd(Gxmedio_paso_t);
	gxmediopasot_menos_mt = _mm256_set1_pd(-Gxmedio_paso_t);
	val_coor_mt = (double*)&coor_despues_mt;
	val_vi_mt = (double*)&vi_despues_mt;
	val_v_mt = (double*)&v_despues_mt;
#endif // AVX2

	// Calcular

	iteracion = 0;
	while (iteracion < hasta) {
		if (!IncrementaTMono()) return;
		if (nparticulas == 1) break;
		iteracion++;
		if (iteracion % frecu_monitorizar == 0) {
			Seguimiento(iteracion);
		}
		if (iteracion == hasta) break;
		if (iteracion % frecu_salidas == 0) {
			EscribeDatos(nsalida++);
		}
		if (CANCELADO) break;
	}
	Terminacion(iteracion, hasta, "Procesamiento mono tarea");
}

bool IncrementaTMono() {

	// Calcular la aceleración de cada partícula

	Aceleraciones();

	// Salto de rana

#ifndef AVX2
	for (int n = 0; n < nparticulas; n++) {
		if (masas[n] == 0) break;
		for (int k = 0; k < 3; k++) {
			velintermedias[k][n] = velintermedias[k][n] + aceleraciones[k][n] * Gxpaso_t;
			coordenadas[k][n] = coordenadas[k][n] + velintermedias[k][n] * paso_t;
			velocidades[k][n] = velintermedias[k][n] - aceleraciones[k][n] * Gxmedio_paso_t;
		}
	}
#else
	for (int n = 0; n < npm4_mt; n += bloque_mt) {
		if (masas[n] == 0) break;
		for (int k = 0; k < 3; k++) {

			/*coor_antes_mt = _mm256_load_pd(&coordenadas[k][n]);
			vi_antes_mt = _mm256_load_pd(&velintermedias[k][n]);
			a_antes_mt = _mm256_load_pd(&aceleraciones[k][n]);*/

			coor_antes_mt = _mm256_setr_pd(coordenadas[k][n], coordenadas[k][n + 1], coordenadas[k][n + 2], coordenadas[k][n + 3]);
			vi_antes_mt = _mm256_setr_pd(velintermedias[k][n], velintermedias[k][n + 1], velintermedias[k][n + 2], velintermedias[k][n + 3]);
			a_antes_mt = _mm256_setr_pd(aceleraciones[k][n], aceleraciones[k][n + 1], aceleraciones[k][n + 2], aceleraciones[k][n + 3]);
			vi_despues_mt = _mm256_fmadd_pd(a_antes_mt, gxpasot_mt, vi_antes_mt);
			coor_despues_mt = _mm256_fmadd_pd(vi_despues_mt, pasot_mt, coor_antes_mt);
			v_despues_mt = _mm256_fmadd_pd(a_antes_mt, gxmediopasot_menos_mt, vi_despues_mt);

			for (int j = 0; j < bloque_mt; j++) {
				velintermedias[k][n + j] = val_vi_mt[j];
				coordenadas[k][n + j] = val_coor_mt[j];
				velocidades[k][n + j] = val_v_mt[j];
			}
		}
	}
	if (npm4_mt < npt_mt) {
		for (int n = npm4_mt; n < npt_mt; n++) {
			if (masas[n] == 0) break;
			for (int k = 0; k < 3; k++) {
				velintermedias[k][n] = velintermedias[k][n] + aceleraciones[k][n] * Gxpaso_t;
				coordenadas[k][n] = coordenadas[k][n] + velintermedias[k][n] * paso_t;
				velocidades[k][n] = velintermedias[k][n] - aceleraciones[k][n] * Gxmedio_paso_t;
			}
		}
	}
#endif // AVX2
	int j;
	double dx;
	double dy;
	double dz;
	for (int i = 0; i < nseguir; i++) {
		j = particulas_seguir[i];
		dx = velintermedias[0][j] * paso_t;
		dy = velintermedias[1][j] * paso_t;
		dz = velintermedias[2][j] * paso_t;
		recorrido[i] += sqrt(dx * dx + dy * dy + dz * dz);
	}
	return true;
}

void PreparaMP(char* ficheromuestra, unsigned long hasta) {

	// Comprobar la relación entre el número de hilos y el de partículas

	int gh = nparticulas / n_hilos;
	if (gh * n_hilos != nparticulas) {
		printf("Cancelado. El numero de partículas debe ser multiplo del numero de hilos");
		return;
	}
	if (!Arranque(ficheromuestra, hasta)) return;

	unsigned long iteracion = ProcesoMP(hasta, gh);

	if (nparticulas == 1) {
		printf("\nEl universo ha coplapsado\n");
	}
	Terminacion(iteracion, hasta, "Procesamiento multi tarea");
}

unsigned long ProcesoMP(unsigned long hasta, int gh) {

	// Antes de crear los hilos hay que calcular las velocidades intermedias y luego poner a cero las aceleraciones porque arrancan según se crean

	// Crear las secciones críticas

	if (!InitializeCriticalSectionAndSpinCount(&seccion_fin_hilos, 0x00000400)) {
		printf("Error [%d] creando la seccion critica fin_hilos\n", GetLastError());
		return 0;
	}
	if (!InitializeCriticalSectionAndSpinCount(&seccion_colision, 0x00000400)) {
		printf("Error [%d] creando la seccion critica fusión\n", GetLastError());
		return 0;
	}
#ifdef HILO_VECINAS
	if (u_vecinas == 1) {
		// Crear el semáforo para la tabla de vecinas abierto y su hilo

		semaforo_vecinas = CreateSemaphore(NULL, 1, 1, NULL);
		if (semaforo_vecinas == NULL)
		{
			printf("Error al crear el semaforo de vecinas: %d\n", GetLastError());
			return 0;
		}
		hilo_vecinas = (HANDLE)_beginthreadex(0, 0, HiloVecinas, NULL, 0, 0);
		Sleep(50);
		if (hilo_vecinas == 0) {
			printf("Error al crear el hilo de vecinas");
			return 0;
		}
	}
#endif // HILO_VECINAS

	// Crear el semáforo para iterar cerrado

	semaforo_itera = CreateSemaphore(NULL, 0, 1, NULL);
	if (semaforo_itera == NULL)
	{
		printf("Error al crear el semaforo de iteracion: %d\n", GetLastError());
		return 0;
	}

	// Crear los semaforos para los hilos abiertos

	hilos_terminados = 0;
	for (int i = 0; i < n_hilos; i++) {
		semaforo_hilo[i] = CreateSemaphore(NULL, 1, 1, NULL);
		if (semaforo_hilo[i] == NULL)
		{
			printf("Error al crear el semaforo para los hilos: %d\n", GetLastError());
			return 0;
		}
	}

	// Arranca los hilos con su semáforo abierto

	int parametros_hilo[MAX_HILOS * 3];
	int conth = 0;
	for (int i = 0; i < n_hilos; i++) {
		parametros_hilo[conth++] = i;
		parametros_hilo[conth++] = i * gh;
		if (i == n_hilos - 1) {
			parametros_hilo[conth++] = nparticulas;
		}
		else {
			parametros_hilo[conth++] = i * gh + gh;
		}
	}
	for (int i = 0; i < n_hilos; i++) {
		hilos[i] = (HANDLE)_beginthreadex(0, 0, HiloAceleraciones, parametros_hilo + 3 * (INT_PTR)i, 0, 0);

		// Para que los hilos comienzen a ejecutarse en orden

		Sleep(50);

		if (hilos[i] == 0) {
			printf("Error al crear un hilo");
			return 0;
		}
	}
	int ns;
	double dx;
	double dy;
	double dz;
#ifdef AVX2
	const int bloque = 4;
	int npt = nparticulas;
	int npm4 = (npt / bloque) * bloque;
	int nm1;
	int nm2;
	int nm3;
	__m256d pasot = _mm256_set1_pd(paso_t);
	__m256d gxpasot = _mm256_set1_pd(Gxpaso_t);
	__m256d gxmediopasot = _mm256_set1_pd(Gxmedio_paso_t);
	__m256d gxmediopasot_menos = _mm256_set1_pd(-Gxmedio_paso_t);
	__m256d coor_antes;
	__m256d vi_antes;
	__m256d a_antes;
	__m256d coor_despues;
	__m256d vi_despues;
	__m256d v_despues;
	double* val_coor_despues = (double*)&coor_despues;
	double* val_vi_despues = (double*)&vi_despues;
	double* val_v_despues = (double*)&v_despues;
#endif // AVX2

	// Arranca el ciclo temporal con el semáforo iterar cerrado. Lo abre el último hilo en terminar

	DWORD espera;
	iteracion = 0;
	while (iteracion < hasta) {
		espera = WaitForSingleObject(semaforo_itera, INFINITE);
		switch (espera)
		{
		case WAIT_OBJECT_0:

			// Los hilos han terminado de calcular las aceleraciones

			// Salto de rana

#ifndef AVX2
			for (int n = 0; n < nparticulas; n++) {
				if (masas[n] == 0) break;
				for (int k = 0; k < 3; k++) {
					velintermedias[k][n] = velintermedias[k][n] + aceleraciones[k][n] * Gxpaso_t;
					coordenadas[k][n] = coordenadas[k][n] + velintermedias[k][n] * paso_t;
					velocidades[k][n] = velintermedias[k][n] - aceleraciones[k][n] * Gxmedio_paso_t;
				}
			}
#else
			for (int n = 0; n < npm4; n += bloque) {
				if (masas[n] == 0) break;
				coor_antes = _mm256_setr_pd(coordenadas[0][n], coordenadas[0][nm1 = n + 1], coordenadas[0][nm2 = n + 2], coordenadas[0][nm3 = n + 3]);
				vi_antes = _mm256_setr_pd(velintermedias[0][n], velintermedias[0][nm1], velintermedias[0][nm2], velintermedias[0][nm3]);
				a_antes = _mm256_setr_pd(aceleraciones[0][n], aceleraciones[0][nm1], aceleraciones[0][nm2], aceleraciones[0][nm3]);
				vi_despues = _mm256_fmadd_pd(a_antes, gxpasot, vi_antes);
				coor_despues = _mm256_fmadd_pd(vi_despues, pasot, coor_antes);
				v_despues = _mm256_fmadd_pd(a_antes, gxmediopasot_menos, vi_despues);
				velintermedias[0][n] = val_vi_despues[0];
				coordenadas[0][n] = val_coor_despues[0];
				velocidades[0][n] = val_v_despues[0];
				velintermedias[0][nm1] = val_vi_despues[1];
				coordenadas[0][nm1] = val_coor_despues[1];
				velocidades[0][nm1] = val_v_despues[1];
				velintermedias[0][nm2] = val_vi_despues[2];
				coordenadas[0][nm2] = val_coor_despues[2];
				velocidades[0][nm2] = val_v_despues[2];
				velintermedias[0][nm3] = val_vi_despues[3];
				coordenadas[0][nm3] = val_coor_despues[3];
				velocidades[0][nm3] = val_v_despues[3];
				coor_antes = _mm256_setr_pd(coordenadas[1][n], coordenadas[1][nm1], coordenadas[1][nm2], coordenadas[1][nm3]);
				vi_antes = _mm256_setr_pd(velintermedias[1][n], velintermedias[1][nm1], velintermedias[1][nm2], velintermedias[1][nm3]);
				a_antes = _mm256_setr_pd(aceleraciones[1][n], aceleraciones[1][nm1], aceleraciones[1][nm2], aceleraciones[1][nm3]);
				vi_despues = _mm256_fmadd_pd(a_antes, gxpasot, vi_antes);
				coor_despues = _mm256_fmadd_pd(vi_despues, pasot, coor_antes);
				v_despues = _mm256_fmadd_pd(a_antes, gxmediopasot_menos, vi_despues);
				velintermedias[1][n] = val_vi_despues[0];
				coordenadas[1][n] = val_coor_despues[0];
				velocidades[1][n] = val_v_despues[0];
				velintermedias[1][nm1] = val_vi_despues[1];
				coordenadas[1][nm1] = val_coor_despues[1];
				velocidades[1][nm1] = val_v_despues[1];
				velintermedias[1][nm2] = val_vi_despues[2];
				coordenadas[1][nm2] = val_coor_despues[2];
				velocidades[1][nm2] = val_v_despues[2];
				velintermedias[1][nm3] = val_vi_despues[3];
				coordenadas[1][nm3] = val_coor_despues[3];
				velocidades[1][nm3] = val_v_despues[3];
				coor_antes = _mm256_setr_pd(coordenadas[2][n], coordenadas[2][nm1], coordenadas[2][nm2], coordenadas[2][nm3]);
				vi_antes = _mm256_setr_pd(velintermedias[2][n], velintermedias[2][nm1], velintermedias[2][nm2], velintermedias[2][nm3]);
				a_antes = _mm256_setr_pd(aceleraciones[2][n], aceleraciones[2][nm1], aceleraciones[2][nm2], aceleraciones[2][nm3]);
				vi_despues = _mm256_fmadd_pd(a_antes, gxpasot, vi_antes);
				coor_despues = _mm256_fmadd_pd(vi_despues, pasot, coor_antes);
				v_despues = _mm256_fmadd_pd(a_antes, gxmediopasot_menos, vi_despues);
				velintermedias[2][n] = val_vi_despues[0];
				coordenadas[2][n] = val_coor_despues[0];
				velocidades[2][n] = val_v_despues[0];
				velintermedias[2][nm1] = val_vi_despues[1];
				coordenadas[2][nm1] = val_coor_despues[1];
				velocidades[2][nm1] = val_v_despues[1];
				velintermedias[2][nm2] = val_vi_despues[2];
				coordenadas[2][nm2] = val_coor_despues[2];
				velocidades[2][nm2] = val_v_despues[2];
				velintermedias[2][nm3] = val_vi_despues[3];
				coordenadas[2][nm3] = val_coor_despues[3];
				velocidades[2][nm3] = val_v_despues[3];
			}
			if (npm4 < npt) {
				for (int n = npm4; n < npt; n++) {
					if (masas[n] == 0) break;
					velintermedias[0][n] = velintermedias[0][n] + aceleraciones[0][n] * Gxpaso_t;
					coordenadas[0][n] = coordenadas[0][n] + velintermedias[0][n] * paso_t;
					velocidades[0][n] = velintermedias[0][n] - aceleraciones[0][n] * Gxmedio_paso_t;
					velintermedias[1][n] = velintermedias[1][n] + aceleraciones[1][n] * Gxpaso_t;
					coordenadas[1][n] = coordenadas[1][n] + velintermedias[1][n] * paso_t;
					velocidades[1][n] = velintermedias[1][n] - aceleraciones[1][n] * Gxmedio_paso_t;
					velintermedias[2][n] = velintermedias[2][n] + aceleraciones[2][n] * Gxpaso_t;
					coordenadas[2][n] = coordenadas[2][n] + velintermedias[2][n] * paso_t;
					velocidades[2][n] = velintermedias[2][n] - aceleraciones[2][n] * Gxmedio_paso_t;
				}
			}
#endif // AVX2
			for (int i = 0; i < nseguir; i++) {
				ns = particulas_seguir[i];
				dx = velintermedias[0][ns] * paso_t;
				dy = velintermedias[1][ns] * paso_t;
				dz = velintermedias[2][ns] * paso_t;
				recorrido[i] += sqrt(dx * dx + dy * dy + dz * dz);
			}
			iteracion++;
			if (iteracion % frecu_monitorizar == 0)
			{
				Seguimiento(iteracion);
			}
			if (iteracion == hasta) break;
			if (iteracion % frecu_salidas == 0) {
				EscribeDatos(nsalida++);
			}
			if (!CANCELADO) {
				if (u_vecinas == 1) {
					if (iteracion % f_vecinas == 0) {
#ifdef HILO_VECINAS
						// Libera el semáforo para la tabla de vecinas

						while (ReleaseSemaphore(semaforo_vecinas, 1, NULL) == FALSE);
#else
						ActualizaTablaVecinas();
#endif // HILO_VECINAS
					}
				}
				if (t_vecinas_enuso != t_vecinas_actualizada) {

					// Cambia la tabla de vecinas más actualizada

					t_vecinas_enuso = t_vecinas_actualizada;
				}

				// Libera los hilos para el cálculo de aceleraciones

				hilos_terminados = 0;
				for (int i = 0; i < n_hilos; i++) {
					while (ReleaseSemaphore(semaforo_hilo[i], 1, NULL) == FALSE);
				}
			}
			break;
		case WAIT_TIMEOUT:
			printf("El hilo %d: supero el tiempo de espera\n", GetCurrentThreadId());
			break;
		}
		if (CANCELADO) break;
	}

	/*

	// Provoca un error en los hilos que siguen activos

	CloseHandle(semaforo_itera);
	for (int i = 0; i < n_hilos; i++) {
		CloseHandle(hilos[i]);
		CloseHandle(semaforo_hilo[i]);
	}
	DeleteCriticalSection(&seccion_fin_hilos);
	DeleteCriticalSection(&seccion_fusion); */

	return iteracion;
}

unsigned int WINAPI HiloAceleraciones(LPVOID param)
{
	int* parametros = (int*)param;
	int indice = parametros[0];
	int desde = parametros[1];
	int hasta = parametros[2];

	// Ciclo perpetuo

	//DWORD resBloqueo;
	DWORD resSemaforo;
	while (TRUE)
	{
		resSemaforo = WaitForSingleObject(semaforo_hilo[indice], INFINITE);
		//if (CANCELADO)return TRUE;
		switch (resSemaforo)
		{
		case WAIT_OBJECT_0:
			AceleracionesTramo(indice, desde, hasta);
			EnterCriticalSection(&seccion_fin_hilos);
			hilos_terminados++;
			if (hilos_terminados == n_hilos) {

				// Libera el semáforo para iterar

				while (ReleaseSemaphore(semaforo_itera, 1, NULL) == FALSE);
				/*if (!ReleaseSemaphore(semaforo_itera, 1, NULL))
				{
					printf("Error [%d] liberando el semáforo para iterar en el hilo: %d: \n", GetLastError(), indice);
					return 0;
				}*/
			}
			LeaveCriticalSection(&seccion_fin_hilos);
			break;
		case WAIT_TIMEOUT:
			printf("El hilo %d supero el tiempo de espera\n", indice);
			break;
		}
	}
	return TRUE;
	return 0;
}

#ifdef HILO_VECINAS
unsigned int WINAPI HiloVecinas(LPVOID param)
{
	// Ciclo perpetuo

	DWORD resSemaforo;
	while (TRUE)
	{
		resSemaforo = WaitForSingleObject(semaforo_vecinas, INFINITE);
		switch (resSemaforo)
		{
		case WAIT_OBJECT_0:
			ActualizaTablaVecinas();
			break;
		case WAIT_TIMEOUT:
			printf("El hilo de vecinas supero el tiempo de espera\n");
			break;
		}
	}
	return TRUE;
	return 0;
}
#endif HILO_VECINAS

void Aceleraciones() {
	double x;
	double y;
	double z;
	double d;
	double d2;
	double f;
	double a;
	for (int i = 0; i < nparticulas; i++) {
		aceleraciones[0][i] = 0.0;
		aceleraciones[1][i] = 0.0;
		aceleraciones[2][i] = 0.0;
	}
	for (int i = 0; i < nparticulas; i++) {
		if (masas[i] == 0.0) break;
		for (int j = i + 1; j < nparticulas; j++) {
			if (masas[j] == 0.0) break;
			x = coordenadas[0][j] - coordenadas[0][i];
			y = coordenadas[1][j] - coordenadas[1][i];
			z = coordenadas[2][j] - coordenadas[2][i];
			d2 = x * x + y * y + z * z;
			if (d2 > radio_corte2) continue;
			d = sqrt(d2);
			if (d < dis_colision) {
				TrataColision(iteracion, i, j, &d2, &d);
			}
			f = d2 * d;
			a = masas[j] / f;
			aceleraciones[0][i] += a * x;
			aceleraciones[1][i] += a * y;
			aceleraciones[2][i] += a * z;
			a = -masas[i] / f;
			aceleraciones[0][j] += a * x;
			aceleraciones[1][j] += a * y;
			aceleraciones[2][j] += a * z;
		}
	}
}

void AceleracionesTramo(int hilo, int desde, int hasta) {
	int nvt;
	int j;
	double x;
	double y;
	double z;
	double d;
	double d2;
	double a;
#ifdef AVX2
	const int bloque = 4;
	int nvmb;
	int n0;
	int n1;
	int n2;
	int n3;
	double aa0;
	double aa1;
	double aa2;
	double aa3;
	__m256d dcxi;
	__m256d dcyi;
	__m256d dczi;
	__m256d dcxj;
	__m256d dcyj;
	__m256d dczj;
	__m256d dcx;
	__m256d dcy;
	__m256d dcz;
	__m256d dc2;
	__m256d dc;
	__m256d dc3;
	__m256d daax;
	__m256d daay;
	__m256d daaz;
	__m256d da;
	__m256d sumah_a;
	double* val_dc = (double*)&dc;
	double* val_dc3 = (double*)&dc3;
	double* val_sumah_a = (double*)&sumah_a;
#endif // AVX2
	for (int i = desde; i < hasta; i++) {

		// Las partículas ficticias están al final de forma consecutiva, terminar al encontrar la primera

		if (masas[i] == 0.0) break;
#ifdef AVX2
		// Poner a cero las aceleraciones sobre la partículas 'i', 'inc' (4) veces ('inc' cumuladores)

		daax = _mm256_setzero_pd();
		daay = _mm256_setzero_pd();
		daaz = _mm256_setzero_pd();

		// Coordenadas de la partícula 'i', 'inc' (4) veces

		dcxi = _mm256_set1_pd(coordenadas[0][i]);
		dcyi = _mm256_set1_pd(coordenadas[1][i]);
		dczi = _mm256_set1_pd(coordenadas[2][i]);

		switch (t_vecinas_enuso)
		{
		case 0:
			nvt = n_vecinas[i][0];
			break;
		case 1:
			nvt = n_vecinas[i][1];
			break;
		default:
			nvt = nparticulas;
			break;
		}
		nvmb = (nvt / bloque) * bloque;
		for (int indj = 0; indj < nvmb; indj += bloque) {
			switch (t_vecinas_enuso)
			{
			case 0:
				n0 = vecinas[i][0][indj];
				n1 = vecinas[i][0][indj + 1];
				n2 = vecinas[i][0][indj + 2];
				n3 = vecinas[i][0][indj + 3];
				break;
			case 1:
				n0 = vecinas[i][1][indj];
				n1 = vecinas[i][1][indj + 1];
				n2 = vecinas[i][1][indj + 2];
				n3 = vecinas[i][1][indj + 3];
				break;
			default:
				n0 = indj;
				n1 = indj + 1;
				n2 = indj + 2;
				n3 = indj + 3;
				break;
			}
			if (masas[n0] == 0.0) break;

			// Calcula distancias de la partícula 'i' respecto a las: 'j' hasta 'j+inc-1'

			dcxj = _mm256_setr_pd(coordenadas[0][n0], coordenadas[0][n1], coordenadas[0][n2], coordenadas[0][n3]);
			dcyj = _mm256_setr_pd(coordenadas[1][n0], coordenadas[1][n1], coordenadas[1][n2], coordenadas[1][n3]);
			dczj = _mm256_setr_pd(coordenadas[2][n0], coordenadas[2][n1], coordenadas[2][n2], coordenadas[2][n3]);

			dcx = _mm256_sub_pd(dcxj, dcxi);
			dcy = _mm256_sub_pd(dcyj, dcyi);
			dcz = _mm256_sub_pd(dczj, dczi);

			dc2 = _mm256_mul_pd(dcx, dcx);
			dc2 = _mm256_fmadd_pd(dcy, dcy, dc2);
			dc2 = _mm256_fmadd_pd(dcz, dcz, dc2);

			dc = _mm256_sqrt_pd(dc2);
			dc3 = _mm256_mul_pd(dc2, dc);

			// Preparar el cálculo del módulo de la aceleración. No se puede vectorizar porque hay que:
			//		Excluir el par 'i'-'i'
			//		Aplicar el radio de corte
			//		Aplicar la distancia de colisión

			if (i == n0 || val_dc[0] > radio_corte) {
				aa0 = 0.0;
			}
			else {
				if (val_dc[0] < dis_colision) {
					TrataColisionHilo(iteracion, i, n0, &d2, &d);
					aa0 = masas[n0] / (d2 * d);
				}
				else {
					aa0 = masas[n0] / val_dc3[0];
				}
			}
			if (i == n1 || masas[n1] == 0.0 || val_dc[1] > radio_corte) {
				aa1 = 0.0;
			}
			else {
				if (val_dc[1] < dis_colision) {
					TrataColisionHilo(iteracion, i, n1, &d2, &d);
					aa1 = masas[n1] / (d2 * d);
				}
				else {
					aa1 = masas[n1] / val_dc3[1];
				}
			}
			if (i == n2 || masas[n2] == 0.0 || val_dc[2] > radio_corte) {
				aa2 = 0.0;
			}
			else {
				if (val_dc[2] < dis_colision) {
					TrataColisionHilo(iteracion, i, n2, &d2, &d);
					aa2 = masas[n2] / (d2 * d);
				}
				else {
					aa2 = masas[n2] / val_dc3[2];
				}
			}
			if (i == n3 || masas[n3] == 0.0 || val_dc[3] > radio_corte) {
				aa3 = 0.0;
			}
			else {
				if (val_dc[3] < dis_colision) {
					TrataColisionHilo(iteracion, i, n3, &d2, &d);
					aa3 = masas[n3] / (d2 * d);
				}
				else {
					aa3 = masas[n3] / val_dc3[3];
				}
			}
			da = _mm256_setr_pd(aa0, aa1, aa2, aa3);

			// Proyecta el módulo de la 'aceleración/r' sobre los ejes, además acumula las componentes de la aceleración
			// resultantes sobre los 'inc' (4) acumuladores de la partícula 'i'

			daax = _mm256_fmadd_pd(da, dcx, daax);
			daay = _mm256_fmadd_pd(da, dcy, daay);
			daaz = _mm256_fmadd_pd(da, dcz, daaz);
		}

		// Suma los 'inc' acumuladores de la partícula 'i'

		sumah_a = _mm256_hadd_pd(daax, daax);
		aceleraciones[0][i] = val_sumah_a[0] + val_sumah_a[2];
		sumah_a = _mm256_hadd_pd(daay, daay);
		aceleraciones[1][i] = val_sumah_a[0] + val_sumah_a[2];
		sumah_a = _mm256_hadd_pd(daaz, daaz);
		aceleraciones[2][i] = val_sumah_a[0] + val_sumah_a[2];

		if (nvmb < nvt) {

			// El resto de partículas 'j', después de dividir por 'bloque', una a una

			for (int indj = nvmb; indj < nvt; indj++) {
				switch (t_vecinas_enuso)
				{
				case 0:
					j = vecinas[i][0][indj];
					break;
				case 1:
					j = vecinas[i][1][indj];
					break;
				default:
					j = indj;
					break;
				}
				if (i == j) continue;
				if (masas[j] == 0.0) break;
				x = coordenadas[0][j] - coordenadas[0][i];
				y = coordenadas[1][j] - coordenadas[1][i];
				z = coordenadas[2][j] - coordenadas[2][i];
				d2 = x * x + y * y + z * z;
				if (d2 > radio_corte2) continue;
				d = sqrt(d2);
				if (d < dis_colision) {
					TrataColisionHilo(iteracion, i, j, &d2, &d);
				}
				a = masas[j] / (d2 * d);
				aceleraciones[0][i] += a * x;
				aceleraciones[1][i] += a * y;
				aceleraciones[2][i] += a * z;
			}
		}
#else
		aceleraciones[0][i] = 0.0;
		aceleraciones[1][i] = 0.0;
		aceleraciones[2][i] = 0.0;
		switch (t_vecinas_enuso)
		{
		case 0:
			nvt = n_vecinas[i][0];
			break;
		case 1:
			nvt = n_vecinas[i][1];
			break;
		default:
			nvt = nparticulas;
			break;
		}
		for (int indj = 0; indj < nvt; indj++) {
			switch (t_vecinas_enuso)
			{
			case 0:
				j = vecinas[i][0][indj];
				break;
			case 1:
				j = vecinas[i][1][indj];
				break;
			default:
				j = indj;
				break;
			}
			if (i == j) continue;
			if (masas[j] == 0.0) break;
			x = coordenadas[0][j] - coordenadas[0][i];
			y = coordenadas[1][j] - coordenadas[1][i];
			z = coordenadas[2][j] - coordenadas[2][i];
			d2 = x * x + y * y + z * z;
			if (d2 > radio_corte2) continue;
			d = sqrt(d2);
			if (d < dis_colision) {
				TrataColision(iteracion, i, j, &d2, &d);
			}
			a = masas[j] / (d2 * d);
			aceleraciones[0][i] += a * x;
			aceleraciones[1][i] += a * y;
			aceleraciones[2][i] += a * z;
		}
#endif // AVX2
	}
}

void ActualizaTablaVecinas() {
	int indice;
	switch (t_vecinas_enuso)
	{
	case 0:
		indice = 1;
		break;
	case 1:
		indice = 0;
		break;
	default:
		indice = 0;
		break;
		}
	double x;
	double y;
	double z;
	double d2;
	int n;
	for (int i = 0; i < nparticulas; i++) {
		n = 0;
		if (masas[i] == 0.0) break;
		for (int j = 0; j < nparticulas; j++) {
			if (i == j) continue;
			if (masas[j] == 0.0) break;
			x = coordenadas[0][i] - coordenadas[0][j];
			y = coordenadas[1][i] - coordenadas[1][j];
			z = coordenadas[2][i] - coordenadas[2][j];
			d2 = x * x + y * y + z * z;
			if (d2 > radio_vecinas2) continue;
			vecinas[i][indice][n++] = j;
			if (n == MAX_VECINAS) {
				printf("Superado el limite de vecinas. Anulado el uso de la tabla.\n");
				fprintf(elog, "\nSuperado el limite de vecinas. Anulado el uso de la tabla.\n\n");
				t_vecinas_actualizada = t_vecinas_enuso = -1;
				u_vecinas = 0;
				return;
			}
		}
		n_vecinas[i][indice] = n;
		if (n > max_absoluto_vecinas) {
			max_absoluto_vecinas = n;
		}
	}

	// Tabla actualizada, se empezará a usar una vez hayan terminado todos los hilos que calculan aceleraciones

	t_vecinas_actualizada = indice;
	}

void TrataColision(unsigned long iteracion, int i, int j, double* pd2, double* pd) {
	
	// 'pd2' y 'pd' son para devolver datos, los valores recibidos no importan

	double x = coordenadas[0][i] - coordenadas[0][j];
	double y = coordenadas[1][i] - coordenadas[1][j];
	double z = coordenadas[2][i] - coordenadas[2][j];
	(*pd2) = x * x + y * y + z * z;
	(*pd) = sqrt(*pd2);
	if ((*pd) >= dis_colision) return;
	double ep;
	double eci;
	double ecj;
	if ((*pd) < dis_critica) {

		// Se modifican las coordenadas para separarlas más allá de la distancia crítica

		ncriticas++;
		fprintf(ecol, "%14lu;%14d;%14d;%14d;%14.7e CRITICA antes\n", iteracion, ncriticas, i, j, *pd);
		MuestraDatosParticula(ecol, i);
		MuestraDatosParticula(ecol, j);
		fflush(ecol);
		double cc[3][2];
		for (int k = 0; k < 3; k++) {
			cc[k][0] = coordenadas[k][i];
			cc[k][1] = coordenadas[k][j];
		}
		double nd2;
		double nd;
		while ((*pd) < dis_critica) {

			// Se desplazan en la dirección que llevaban

			for (int k = 0; k < 3; k++) {
				cc[k][0] += velintermedias[k][i] * paso_t;
				cc[k][1] += velintermedias[k][j] * paso_t;
			}
			x = cc[0][0] - cc[0][1];
			y = cc[1][0] - cc[1][1];
			z = cc[2][0] - cc[2][1];
			nd2 = x * x + y * y + z * z;
			nd = sqrt(nd2);
			if (abs(nd - *pd) < dis_critica / 10.0) {

				// Forzar la separación. Probablemente se mueven en la misma dirección.

				for (int k = 0; k < 3; k++) {
					cc[k][0] += dis_critica / 2.0;
					cc[k][1] -= dis_critica / 2.0;
				}
				x = cc[0][0] - cc[0][1];
				y = cc[1][0] - cc[1][1];
				z = cc[2][0] - cc[2][1];
				nd2 = x * x + y * y + z * z;
				nd = sqrt(nd2);
			}
			(*pd) = nd;
		}
		(*pd2) = nd2;
		for (int k = 0; k < 3; k++) {
			coordenadas[k][i] = cc[k][0];
			coordenadas[k][j] = cc[k][1];
		}
		ep = -G * masas[i] * masas[j] / *pd;
		eci = masas[i] * (velintermedias[0][i] * velintermedias[0][i] + velintermedias[1][i] * velintermedias[1][i] + velintermedias[2][i] * velintermedias[2][i]) / 2.0;
		ecj = masas[j] * (velintermedias[0][j] * velintermedias[0][j] + velintermedias[1][j] * velintermedias[1][j] + velintermedias[2][j] * velintermedias[2][j]) / 2.0;
		fprintf(ecol, "%14lu;%14d;%14d;%14d;%14.7e;%14.7e;%14.7e;%14.7e CRITICA despues\n", iteracion, ncriticas, i, j, *pd, eci, ecj, ep);
		MuestraDatosParticula(ecol, i);
		MuestraDatosParticula(ecol, j);
		fflush(ecol);
	}

	// Si se están aproximando (disminuye la distancia entre ellas) se intercambiam
	// El problema con un segundo hilo no se presenta porque una vez intercambiadas se alejan

	x += (velintermedias[0][i] - velintermedias[0][j]) * paso_t;
	y += (velintermedias[1][i] - velintermedias[1][j]) * paso_t;
	z += (velintermedias[2][i] - velintermedias[2][j]) * paso_t;
	double df = sqrt(x * x + y * y + z * z);
	if (df > (*pd)) {

		// Se están alejando

		return;
	}

	// Se intercambian las posiciones pero se mantienen las velocidades.
	// antes        a -->  <-- b  
	// después  <-- b          a -->

	ncolisiones++;
	double val;
	for (int k = 0; k < 3; k++) {
		val = coordenadas[k][i];
		coordenadas[k][i] = coordenadas[k][j];
		coordenadas[k][j] = val;
	}
	ep = -G * masas[i] * masas[j] / *pd;
	eci = masas[i] * (velintermedias[0][i] * velintermedias[0][i] + velintermedias[1][i] * velintermedias[1][i] + velintermedias[2][i] * velintermedias[2][i]) / 2.0;
	ecj = masas[j] * (velintermedias[0][j] * velintermedias[0][j] + velintermedias[1][j] * velintermedias[1][j] + velintermedias[2][j] * velintermedias[2][j]) / 2.0;
	fprintf(ecol, "%14lu;%14d;%14d;%14d;%14.7e;%14.7e;%14.7e;%14.7e;%14.7e\n", iteracion, ncolisiones, i, j, *pd, df, eci, ecj, ep);
	MuestraDatosParticula(ecol, i);
	MuestraDatosParticula(ecol, j);
	fflush(ecol);
}

void TrataColisionHilo(unsigned long iteracion, int i, int j, double* pd2, double* pd) {
	EnterCriticalSection(&seccion_colision);
	TrataColision(iteracion, i, j, pd2, pd);
	LeaveCriticalSection(&seccion_colision);
}

void LiberaMemoria() {
	t_vecinas_actualizada = t_vecinas_enuso = -1;
	max_absoluto_vecinas = 0;
	if (u_vecinas == 1 && vecinas != NULL) {
		free(vecinas);
	}
#ifdef MEMORIA_DINAMICA

	// PENDIENTE DE REVISAR !!!!

	if (masas != NULL) free(masas);
	if (coordenadas != NULL) {
		for (int i = 0; i < 3; i++) free(coordenadas[i]);
		free(coordenadas);
	}
	if (velocidades != NULL) {
		for (int i = 0; i < 3; i++) free(velocidades[i]);
		free(velocidades);
	}
	if (velintermedias != NULL) {
		for (int i = 0; i < 3; i++) free(velintermedias[i]);
		free(velintermedias);
	}
	if (aceleraciones != NULL) {
		for (int i = 0; i < 3; i++) free(aceleraciones[i]);
		free(aceleraciones);
	}
	if (u_vecinas == 1) {
		if (n_vecinas != NULL) {
			for (int i = 0; i < nparticulas; i++) {
				free(n_vecinas[i]);
			}
			free(n_vecinas);
		}
	}
#endif // MEMORIA_DINAMICA
}

bool AsignaMemoria(int nparticulas) {
	t_vecinas_actualizada = t_vecinas_enuso = -1;
	max_absoluto_vecinas = 0;
	if (u_vecinas == 1) {
		size_t ind1 = nparticulas;
		size_t ind2 = 2;
		size_t ind3 = MAX_VECINAS;

		/* memoria =
			  ind1                punteros 'int**' para [i]
		   +  ind1 * ind2         punteros 'int*'  para [i][j]
		   +  ind1 * ind2 * ind3  datos 'int'
		*/

		size_t memoria = ind1 * sizeof(int**) + ind1 * ind2 * sizeof(int*) + ind1 * ind2 * ind3 * sizeof(int);
		vecinas = (int***)malloc(memoria);
		if (vecinas) {
			for (size_t i = 0; i < ind1; ++i) {
				if (i >= memoria - sizeof(int**)) return false; // Es innecesario, es para evitar el aviso 6386 del "code analyzer'
				vecinas[i] = (int**)(vecinas + ind1) + i * ind2;
				for (size_t j = 0; j < ind2; ++j) {
					vecinas[i][j] = (int*)(vecinas + ind1 + ind1 * ind2) + i * ind2 * ind3 + j * ind3;
				}
			}
		}
}
#ifdef MEMORIA_DINAMICA

	// PENDIENTE DE REVISAR !!!!

	masas = (double*)malloc(nparticulas * sizeof(double));
	size_t t_fila = nparticulas * sizeof(double);
	size_t t_3filas = 3 * t_fila;
	coordenadas = (double**)malloc(t_3filas);
	velocidades = (double**)malloc(t_3filas);
	velintermedias = (double**)malloc(t_3filas);
	aceleraciones = (double**)malloc(t_3filas);
	for (int i = 0; i < 3; i++) {
		coordenadas[i] = (double*)malloc(t_fila);
		velocidades[i] = (double*)malloc(t_fila);
		velintermedias[i] = (double*)malloc(t_fila);
		aceleraciones[i] = (double*)malloc(t_fila);
	}
	if (u_vecinas == 1) {
		// Espacio para 'nparticulas' punteros de tipo 'int*'

		n_vecinas = (int**)malloc(nparticulas * sizeof(int*));

		for (int i = 0; i < nparticulas; i++) {
			n_vecinas[i] = (int*)malloc(MAX_VECINAS * sizeof(int));
		}
	}
#endif // MEMORIA_DINAMICA
	return true;
}

/*
bool FusionaMono(unsigned long iteracion, int i, int j) {
	if (masas[i] == 0.0 || masas[j] == 0.0) return false;

	ncolisiones++;
	double v1;
	double v2;

	double eci = Ecinetica();
	double epri;
	double ndentro;
	double epi = Epotencial(&epri, &ndentro);
	double vt;
	double sm = masas[i] + masas[j];
	double x = coordenadas[0][i] - coordenadas[0][j];
	double y = coordenadas[1][i] - coordenadas[1][j];
	double z = coordenadas[2][i] - coordenadas[2][j];
	double d = sqrt(x * x + y * y + z * z);
	double ep = -G * masas[i] * masas[j] / d;

	fprintf(ecol, "%14llu;%14d;%14d;%14d;%14.7e;%14.7e\n", iteracion, ncolisiones, i, j, d, ep);
	MuestraDatosParticula(ecol, i);
	MuestraDatosParticula(ecol, j);
	fflush(elog);

	// Velocidad del conjunto de las dos partículas conservado la energía cinética de ambas

	double ec = 0.0;
	for (int k = 0; k < 3; k++) {
		v1 = velintermedias[k][i];
		v2 = velintermedias[k][j];
		vt = masas[i] * v1 * v1 + masas[j] * v2 * v2;
		velintermedias[k][i] = sqrt(vt / sm);

		v1 = velocidades[k][i];
		v2 = velocidades[k][j];
		vt = masas[i] * v1 * v1 + masas[j] * v2 * v2;
		velocidades[k][i] = sqrt(vt / sm);

		ec += vt;
	}
	ec /= 2.0;

	double uv = 0.0;
	for (int k = 0; k < 3; k++) {
		uv += velocidades[k][i] * velocidades[k][i];
	}
	uv = sqrt(uv);

	// Energia potencial del sistema antes de agrupar

	double eprs;
	double eps = Epotencial(&eprs, &ndentro);

	// Agrupar la masa y las coordenadas de 'i' y 'j'

	masas[i] = sm;
	for (int k = 0; k < 3; k++) {
		coordenadas[k][i] = (coordenadas[k][i] + coordenadas[k][j]) / 2;
	}

	// desactivar 'j'

	masas[j] = false;

	// Nueva energía potencial del sistema

	double nepr;
	double neps = Epotencial(&nepr, &ndentro);

	// La diferencia se debe transformar en energía cinética.
	// Calculamos la nueba velocidad que corresponde a:
	//		la suma de la energía cinética de ambas partículas
	//		más la variacion de ep del sistema

	double nueva_v;
	if (ec + neps - eps < 0.0) {

		// No es posible conservar la energía, lo mas que puede hacer la partícula fusionada es quedarse quieta

		nueva_v = 0.0;
	}
	else {
		nueva_v = sqrt(2.0 * (ec + neps - eps) / sm);
	}
	double uc;
	double nec = 0.0;
	for (int k = 0; k < 3; k++) {
		uc = velocidades[k][i] / uv;
		velintermedias[k][i] = velocidades[k][i] = nueva_v * uc;
		nec += velocidades[k][i] * velocidades[k][i];
	}
	nec *= (sm / 2.0);
	fprintf(ecol, "Var Ep Sist. %14.7e %14.7e %14.7e\n", eps, neps, nueva_v);
	fprintf(ecol, "Par ini      %14.7e %14.7e\n", ec, ep);
	fprintf(ecol, "Par fin      %14.7e %14.7e\n", nec, 0.0);
	ec = Ecinetica();
	double epr;
	ep = Epotencial(&epr, &ndentro);
	fprintf(ecol, "Sistema ini  %14.7e %14.7e %14.7e\n", eci, epi, eci + epi);
	fprintf(ecol, "Sistema fin  %14.7e %14.7e %14.7e\n", ec, ep, ec + ep);
	fflush(ecol);
	return true;
}

bool Fusiona(unsigned long iteracion, int i, int j) {
	if (masas[i] == 0.0 || masas[j] == 0.0) return false;
	EnterCriticalSection(&seccion_colision);
	ncolisiones++;
	double v1;
	double v2;

	double eci = Ecinetica();
	double epri;
	double ndentro;
	double epi = Epotencial(&epri, &ndentro);
	double vt;
	double sm = masas[i] + masas[j];
	double x = coordenadas[0][i] - coordenadas[0][j];
	double y = coordenadas[1][i] - coordenadas[1][j];
	double z = coordenadas[2][i] - coordenadas[2][j];
	double d = sqrt(x * x + y * y + z * z);
	double ep = -G * masas[i] * masas[j] / d;

	fprintf(ecol, "%14llu;%14d;%14d;%14d;%14.7e;%14.7e\n", iteracion, ncolisiones, i, j, d, ep);
	MuestraDatosParticula(ecol, i);
	MuestraDatosParticula(ecol, j);
	fflush(elog);

	// Velocidad del conjunto de las dos partículas conservado la energía cinética de ambas

	double ec = 0.0;
	for (int k = 0; k < 3; k++) {
		v1 = velintermedias[k][i];
		v2 = velintermedias[k][j];
		vt = masas[i] * v1 * v1 + masas[j] * v2 * v2;
		velintermedias[k][i] = sqrt(vt / sm);

		v1 = velocidades[k][i];
		v2 = velocidades[k][j];
		vt = masas[i] * v1 * v1 + masas[j] * v2 * v2;
		velocidades[k][i] = sqrt(vt / sm);

		ec += vt;
	}
	ec /= 2.0;

	double uv = 0.0;
	for (int k = 0; k < 3; k++) {
		uv += velocidades[k][i] * velocidades[k][i];
	}
	uv = sqrt(uv);

	// Energia potencial del sistema antes de agrupar

	double eprs;
	double eps = Epotencial(&eprs, &ndentro);

	// Agrupar la masa y las coordenadas de 'i' y 'j'

	masas[i] = sm;
	for (int k = 0; k < 3; k++) {
		coordenadas[k][i] = (coordenadas[k][i] + coordenadas[k][j]) / 2;
	}

	// desactivar 'j'

	masas[j] = false;

	// Nueva energía potencial del sistema

	double neprs;
	double neps = Epotencial(&neprs, &ndentro);

	// La diferencia se debe transformar en energía cinética.
	// Calculamos la nueba velocidad que corresponde a:
	//		la suma de la energía cinética de ambas partículas
	//		más la variacion de ep del sistema

	double nueva_v;
	if (ec + neps - eps < 0.0) {

		// No es posible conservar la energía, lo mas que puede hacer la partícula fusionada es quedarse quieta

		nueva_v = 0.0;
	}
	else {
		nueva_v = sqrt(2.0 * (ec + neps - eps) / sm);
	}
	double uc;
	double nec = 0.0;
	for (int k = 0; k < 3; k++) {
		uc = velocidades[k][i] / uv;
		velintermedias[k][i] = velocidades[k][i] = nueva_v * uc;
		nec += velocidades[k][i] * velocidades[k][i];
	}
	nec *= (sm / 2.0);
	fprintf(ecol, "Var Ep Sist. %14.7e %14.7e %14.7e\n", eps, neps, nueva_v);
	fprintf(ecol, "Par ini      %14.7e %14.7e\n", ec, ep);
	fprintf(ecol, "Par fin      %14.7e %14.7e\n", nec, 0.0);
	ec = Ecinetica();
	double epr;
	ep = Epotencial(&epr, &ndentro);
	fprintf(ecol, "Sistema ini  %14.7e %14.7e %14.7e\n", eci, epi, eci + epi);
	fprintf(ecol, "Sistema fin  %14.7e %14.7e %14.7e\n", ec, ep, ec + ep);
	fflush(ecol);
	LeaveCriticalSection(&seccion_colision);
	return true;
}
*/