
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#define _USE_MATH_DEFINES
#define M_PI  3.14159265358979323846  /* pi */
void memcpy_image(int** src, int** dst, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			dst[i][j] = src[i][j];
		}
	}
}
class DCT {
private:
	int N = 8;//Block Size
	int QP = 40;
	int Cols;// y:1280,  UV:650
	int Rows;//Y: 640 ,UV: 360
	int **Pixel;//[cols][rows],whole image pixel, Y[1280][720], UV[640][360]
	int **input_Pixel;
	int **transfrom_Pixel;
	int **Reconstruct;
	double C_zero;
	double C_other;
public:
	DCT(unsigned char* &, int,int,double,double);
	void make_DCT();
	void process_DCT(bool);
	int get_DCT(int,int, int, int);
	int get_RDCT(int, int, int, int);
	void MSE();
	void save_quant(std::string output);
};
int DCT::get_DCT(int y, int x, int dy, int dx) {
	double C_y = (dy == 0 ? C_zero : C_other);
	double C_x = (dx == 0 ? C_zero : C_other);
	double res = 0;
	int target_x = x + dx;
	int target_y = y + dy;
	for (int bx = 0; bx < N; bx++) {
		for (int by = 0; by < N; by++) {
			int ref_y = y + by;
			int ref_x = x + bx;
			res += (double)input_Pixel[ref_y][ref_x] * cos((double)((2 * (double)bx+ 1.0)*dx*M_PI) / (2 * N)) *  cos(((2 * (double)by+1.0)*(double)dy*M_PI) / (2 * N));
		}
	}
	return (int)((double)(res * C_y * C_x)/QP);
}
int DCT::get_RDCT(int y, int x, int dy, int dx) {
	double res = 0;
	int target_x = x + dx;
	int target_y = y + dy;
	for (int bx = 0; bx < N; bx++) {
		for (int by = 0; by < N; by++) {
			int ref_y = y + by;
			int ref_x = x + bx;
			double C_y = (by == 0 ? C_zero : C_other);
			double C_x = (bx == 0 ? C_zero : C_other);
			res += ((double)input_Pixel[ref_y][ref_x] * (double)QP * cos(((2 * (double)dx+ 1.0)*(double)bx*M_PI) / (2 * N)) *  cos(((double)(2 * dy+ 1.0)*(double)by*M_PI) / (2 * N))) * C_x * C_y;
		}
	}
	return (int)res;
}
void DCT::process_DCT(bool is_inverse) {
	for (int i = 0; i < Rows; i+=8) {
		for (int j = 0; j < Cols; j+=8) {
			// PROCESS
			for (int dx = 0; dx < N; dx++) {
				for (int dy = 0; dy < N; dy++) {
					if(!is_inverse) transfrom_Pixel[i + dy][j + dx] = get_DCT(i,j,dy,dx);
					else transfrom_Pixel[i + dy][j + dx] = get_RDCT(i, j, dy, dx);
				}
			}
		}
	}
}
void DCT::MSE() {
	double sum = 0;
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			sum += (Reconstruct[i][j] - Pixel[i][j])*(Reconstruct[i][j] - Pixel[i][j]);
		}
	}
	printf("%lf ", sum/((double)Cols*Rows));
}
void DCT::make_DCT() {
	memcpy_image(Pixel, input_Pixel, Rows, Cols);
	process_DCT(false);
	memcpy_image(transfrom_Pixel, input_Pixel, Rows, Cols);
	process_DCT(true);
	memcpy_image(transfrom_Pixel, Reconstruct, Rows, Cols);
	MSE();
}
void DCT::save_quant(std::string s) {
	memcpy_image(Pixel, input_Pixel, Rows, Cols);
	std::string out_name = "QT_" + std::to_string(QP) + "_" + s;
	std::ofstream output(out_name, std::ios::out | std::ios::binary);
	process_DCT(false);
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			if (transfrom_Pixel[i][j] > 255) transfrom_Pixel[i][j] = 255;
			unsigned char c = transfrom_Pixel[i][j];
			output.write((char*)&c, 1);
		}
	}
	output.close();
}
DCT::DCT(unsigned char* &original_input, int Cols, int Rows, double C_zero, double C_other) {
	this->Cols = Cols;
	this->Rows = Rows;
	this->C_zero = C_zero;
	this->C_other = C_other;
	Pixel = new int*[Rows];
	for (int i = 0; i < Rows; i++) Pixel[i] = new int[Cols];
	Reconstruct = new int*[Rows];
	for (int i = 0; i < Rows; i++) Reconstruct[i] = new int[Cols];
	input_Pixel = new int*[Rows];
	for (int i = 0; i < Rows; i++) input_Pixel[i] = new int[Cols];
	transfrom_Pixel = new int*[Rows];
	for (int i = 0; i < Rows; i++) transfrom_Pixel[i] = new int[Cols];

	int idx = 0;
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			Pixel[i][j] = original_input[idx++];
		}
	}
	
	original_input = original_input + idx;
}
int main()
{
	FILE *pFile_Original = NULL;// 오리지널 파일구조체 정의

	pFile_Original = fopen("PeopleOnStreet_1280x720_30_Original.yuv", "rb"); // yuv파일 로드

	int resolution_size = 1280 * 720 * 3 / 2; // 한프레임의 전체 (YUV)를 포함한 해상도 크기
	int y_resolution_size = 1280 * 720; // 한프레임에서 Y값을 가진 
	int uv_resolution_size = 1280 * 720 / 4;

	unsigned char *read_data_Original = new unsigned char[resolution_size]; // Original.yuv 파일의 정보를 담고 있을 배열 공간 정의
	unsigned char *read_data_DCT = new unsigned char[resolution_size]; //Recons.yuv 파일의 정보를 담고 있을 배열 공간의 정의
	int** Y, **U, **V;

	double sum = 0; // Original과 Recon 파일간의 차이를 제곱해 합산하는 변수
	double YMSE; // 

	double p_sum = 0; //  PSNR 값의 합
	double PSNR; // PSNR 값 저장
	size_t n_size_Original;
	size_t n_size_Recons;
	//Original영상파일의 데이터를 읽는 ifelse
	if (pFile_Original == NULL)
	{
		//파일을 읽지 못하는 경우, 에러처리
		fputs("File error", stderr);
		exit(1);
	}
	else
	{
		n_size_Original = fread(read_data_Original, sizeof(unsigned char), resolution_size, pFile_Original);
	}
	DCT Y_Original = DCT(read_data_Original, 1280, 720, sqrt((double)1/8), sqrt((double)2/8));
	//Y_Original.make_DCT();
	DCT U_Original = DCT(read_data_Original, 640, 320, sqrt((double)1 / 8), sqrt((double)2 / 8));
	DCT V_Original = DCT(read_data_Original, 640, 320, sqrt((double)1 / 8), sqrt((double)2 / 8));
	U_Original.save_quant("_data.u");
	V_Original.save_quant("_data.v");

	fclose(pFile_Original); //파일 닫기

	return 0;

}