#include "Descrete_Fourier_Transform.h"

int main()
{
	//мерность векторного пространства
	int N = 512;

	//Z - данные, DFT_Data - прямое DFT, IDFT_Data - обратное DFT
	std::vector<std::complex<double>> Z(N), FFT_Data(N), IFFT_Data(N), DFT_Data(N), IDFT_Data(N);

	//заполнение данных (вещественная гармоника с частотой 100)
	for (int i = 0; i < N; i++)
	{
		Z[i]._Val[0] = cos(2 * PI * i * 500 / N);// +0.01 * cos(2 * PI * i * 510 / N);
		Z[i]._Val[1] = 0;
	}

	//прямое преобразование
	DFT_Data = DFT(Z);
	//быстрое преобразование
	FFT_Data = FFT(Z);

	//обратное преобразование
	IDFT_Data = IDFT(FFT_Data);
	//быстрое обратное преобразование;
	IFFT_Data = IFFT(FFT_Data);

	//вывод результата
	int SETW = 22;

	std::cout << std::left << std::setw(SETW) << "Number"
		<< std::left << std::setw(SETW) << "Re(Z)"
		<< std::left << std::setw(SETW) << "Im(Z)"
		<< std::left << std::setw(SETW) << "Re(FFT_Data)"
		<< std::left << std::setw(SETW) << "Im(FFT_Data)"
		<< std::left << std::setw(SETW) << "Re(DFT_Data)"
		<< std::left << std::setw(SETW) << "Im(DFT_Data)"
		<< std::left << std::setw(SETW) << "Re(IFFT_Data)"
		<< std::left << std::setw(SETW) << "Im(IFFT_Data)"
		<< std::left << std::setw(SETW) << "Re(IDFT_Data)"
		<< std::left << std::setw(SETW) << "Im(IDFT_Data)" << std::endl;

	for (int i = 0; i < 50; i++)
	{
		std::cout << std::left << std::setw(SETW) << i
			<< std::left << std::setw(SETW) << Z[i].real()
			<< std::left << std::setw(SETW) << Z[i].imag()
			<< std::left << std::setw(SETW) << FFT_Data[i].real()
			<< std::left << std::setw(SETW) << FFT_Data[i].imag()
			<< std::left << std::setw(SETW) << DFT_Data[i].real()
			<< std::left << std::setw(SETW) << DFT_Data[i].imag()
			<< std::left << std::setw(SETW) << IFFT_Data[i].real()
			<< std::left << std::setw(SETW) << IFFT_Data[i].imag()
			<< std::left << std::setw(SETW) << IDFT_Data[i].real()
			<< std::left << std::setw(SETW) << IDFT_Data[i].imag() << std::endl;
	}
}