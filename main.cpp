#include "Descrete_Fourier_Transform.h"
#include <chrono>
#include <fstream>

bool equal(const double& x, const double& y)
{
	if (abs(x - y) <= 0.000000001) {
		return true;
	}
	else return false;
}

void Mod_Arg(const std::vector<std::complex<double>>& Data, std::vector<double>& Mod, std::vector<double>& Arg)
{
	for (int i = 0; i < Data.size(); ++i) {
		Mod[i] = sqrt(Data[i]._Val[0] * Data[i]._Val[0] + Data[i]._Val[1] * Data[i]._Val[1]);

		if (!equal(Data[i]._Val[0], 0) && Data[i]._Val[0] > 0) {
			Arg[i] = atan(abs(Data[i]._Val[1] / Data[i]._Val[0]));
		}
		else if (!equal(Data[i]._Val[0], 0) && Data[i]._Val[0] < 0 && Data[i]._Val[1] >= 0) {
			Arg[i] = PI - atan(abs(Data[i]._Val[1] / Data[i]._Val[0]));
		}
		else if (!equal(Data[i]._Val[0], 0) && Data[i]._Val[0] < 0 && !equal(Data[i]._Val[1], 0) && Data[i]._Val[1] < 0) {
			Arg[i] = PI + atan(abs(Data[i]._Val[1] / Data[i]._Val[0]));
		}
		else if (equal(Data[i]._Val[0], 0) && !equal(Data[i]._Val[1], 0) && Data[i]._Val[1] > 0) {
			Arg[i] = PI / 2;
		}
		else if (equal(Data[i]._Val[0], 0) && !equal(Data[i]._Val[1], 0) && Data[i]._Val[1] < 0) {
			Arg[i] = -PI / 2;
		}
	}
}

int main()
{
	//мерность векторного пространства
	int N = 512;

	//Z - данные, DFT_Data - прямое DFT, IDFT_Data - обратное DFT
	std::vector<std::complex<double>> Z(N), Z_n(N), FFT_Data(N), IFFT_Data(N), DFT_Data(N), IDFT_Data(N);
	std::vector<double> Arg(N), Mod(N);

	//заполнение данных 
	for (int i = 0; i < N; i++){
		Z[i]._Val[0] = 23 - 0.77 * cos((2 * PI * i * 231 / N) + PI / 6);
		Z[i]._Val[1] = 0;
		Z_n[i]._Val[0] = cos(2 * PI * i / N) + 0.01 * cos(2 * PI * i * 231 / N);
		Z_n[i]._Val[1] = 0;
	}

	//прямое преобразование
	auto t1 = std::chrono::high_resolution_clock::now();
	DFT_Data = DFT(Z);
	auto t2 = std::chrono::high_resolution_clock::now();
	double time1 = std::chrono::duration<double>(t2 - t1).count();

	//быстрое преобразование
	t1 = std::chrono::high_resolution_clock::now();
	FFT_Data = FFT(Z);
	t2 = std::chrono::high_resolution_clock::now();
	double time2 = std::chrono::duration<double>(t2 - t1).count();

	//обратное преобразование
	t1 = std::chrono::high_resolution_clock::now();
	IDFT_Data = IDFT(FFT_Data);
	t2 = std::chrono::high_resolution_clock::now();
	double time3 = std::chrono::duration<double>(t2 - t1).count();

	//быстрое обратное преобразование;
	t1 = std::chrono::high_resolution_clock::now();
	IFFT_Data = IFFT(FFT_Data);
	t2 = std::chrono::high_resolution_clock::now();
	double time4 = std::chrono::duration<double>(t2 - t1).count();

	std::cout << "DFT Time: " << time1 << std::endl;
	std::cout << "FFT Time: " << time2 << std::endl;
	std::cout << "IDFT Time: " << time3 << std::endl;
	std::cout << "IFFT Time: " << time4 << std::endl;
	std::cout << std::endl;

	Mod_Arg(FFT_Data, Mod, Arg);

	//вывод результата
	int SETW = 22;
	std::cout << "Signal: z(j) = 23 - 0.77*cos(2*pi*j*231/512 + pi/6)" << std::endl;
	std::cout << std::left << std::setw(SETW) << "Number"
	<< std::left << std::setw(SETW) << "Re(Z)"
	<< std::left << std::setw(SETW) << "Re(Z^)"
	<< std::left << std::setw(SETW) << "Im(Z^)"
	<< std::left << std::setw(SETW) << "|Z^|"
	<< std::left << std::setw(SETW) << "arg(Z^)" << std::endl;

	for (int i = 0; i < N; i++){
		if(!equal(Mod[i], 0) || !equal(Arg[i], 0)) {
			std::cout << std::left << std::setw(SETW) << i
			<< std::left << std::setw(SETW) << Z[i].real()
			<< std::left << std::setw(SETW) << FFT_Data[i].real()
			<< std::left << std::setw(SETW) << FFT_Data[i].imag()
			<< std::left << std::setw(SETW) << Mod[i]
			<< std::left << std::setw(SETW) << Arg[i] << std::endl;
		}
	}
	std::cout << std::endl;

	//для задачи 2
	DFT_Data = DFT(Z_n);
	Mod_Arg(DFT_Data, Mod, Arg);

	std::cout << "Signal: z(j) = cos(2*pi*j/N) + 0.01*cos(2*pi*j*231/512)" << std::endl;
	std::cout << std::left << std::setw(SETW) << "Number"
		<< std::left << std::setw(SETW) << "Re(Z)"
		<< std::left << std::setw(SETW) << "Re(Z^)"
		<< std::left << std::setw(SETW) << "Im(Z^)"
		<< std::left << std::setw(SETW) << "|Z^|"
		<< std::left << std::setw(SETW) << "arg(Z^)" << std::endl;

	for (int i = 0; i < N; i++) {
		if (!equal(Mod[i], 0) || !equal(Arg[i], 0)) {
			std::cout << std::left << std::setw(SETW) << i
				<< std::left << std::setw(SETW) << Z[i].real()
				<< std::left << std::setw(SETW) << DFT_Data[i].real()
				<< std::left << std::setw(SETW) << DFT_Data[i].imag()
				<< std::left << std::setw(SETW) << Mod[i]
				<< std::left << std::setw(SETW) << Arg[i] << std::endl;
		}
	}

	DFT_Data[231]._Val[0] = 0;
	DFT_Data[231]._Val[1] = 0;
	DFT_Data[281]._Val[0] = 0;
	DFT_Data[281]._Val[1] = 0;

	IDFT_Data = IDFT(DFT_Data);
	std::ofstream file1("signal_before.txt");
	std::ofstream file2("signal_after.txt");
	for (int i = 0; i < N; ++i) {
		file1 << Z_n[i]._Val[0] << "\n";
		file2 << IDFT_Data[i]._Val[0] << "\n";
		//std::cout << Z_n[i] << " " << IDFT_Data[i] << std::endl;
	}
	file1.close();
	file2.close();
}