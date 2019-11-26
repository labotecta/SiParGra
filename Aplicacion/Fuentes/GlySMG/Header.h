#pragma once
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <signal.h>
#include <string.h>
#include <windows.h>
#include <conio.h>
#include <strsafe.h>
#include <iostream>
#include <tchar.h>
#include <process.h>
#include <psapi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>
#include <limits.h>

#include <thread>
#include <chrono>
#include <random>
#include <algorithm>
#include <direct.h>

#define _USE_MATH_DEFINES // constantes para C
#include <math.h>

// Para detectar las capacidades del procesador

#ifdef _WIN32
//  Windows
#define cpuid(info, x)    __cpuidex(info, x, 0)
#else
//  GCC Intrinsics
#include <cpuid.h>
void cpuid(int info[4], int InfoType) {
	__cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}
#endif

#define AVX2
//#undef AVX2
#ifdef AVX2
#include <immintrin.h>
#endif // AVX2

#define MEMORIA_DINAMICA
#undef MEMORIA_DINAMICA

#define UNIDADES_SMI
#undef UNIDADES_SMI

#define HILO_VECINAS
//#undef HILO_VECINAS
