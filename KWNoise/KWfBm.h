#pragma once
#include "KWNoise.h"

//fBmでノイズを詳細化したテクスチャとして扱うクラス
//template引数 T : ノイズに使うクラス（KWNoiseの派生クラス）
//template引数 N : 内部ノイズに使う解像度（デフォルト128） 
template <class T, int N = 128>
class KWfBm
{
public:
	//初期化
	// dimension (int) : テクスチャの次元
	// resolutions (int*) : テクスチャの解像度．dimension次元で指定
	// time_scale (int) : 時間方向の解像度
	// delta_time (double) : next()のコールごとに進める時間幅
	// octaves (double) : オクターブ数．重ね併せる周波数の数
	// lacunarity (double) : 空隙性[0,1]．重ね併せる周波数の幅
	// hurst_index (double) : ハースト定数[0,1]．フラクタル増分のパラメータ
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	KWfBm(int dimension, int* resolutions, int time_scale, double delta_time,
		double octaves, double lacunarity, double hurst_index, int seed = -1) :
		dim(dimension),
		ts(time_scale),
		dt(delta_time),
		oct(octaves),
		lac(lacunarity),
		H(hurst_index),
		time(0),
		INTERNAL_NOISE_RESOLUTION(N)
	{
		res = new int[dim];
		size = 1;
		for (int i = 0; i < dim; i++)
			size *= (res[i] = resolutions[i]);

		init(seed);
		next();
	}
	~KWfBm() {
		delete res;
		delete aexp;
		delete texture;
		delete noise;
	}

	//ノイズの時間を進める
	void next() {
		double* vec = new double[dim + 1];

		double mx = DBL_MIN;
		double mn = DBL_MAX;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				int step = 1;
				for (int k = 0; k < j; k++)
					step *= res[k];

				vec[j] = ((i / step) % res[j]) / (double)res[j];
			}
			vec[dim] = time * dt / (double)ts;

			double v = texture[i] = get(vec);
			if (mx < v) mx = v;
			if (mn > v) mn = v;
		}

		for (int i = 0; i < size; i++)
		{
			texture[i] = texture[i] / mx;
		}
		
		++time;
		if (time > ts) time = 0;

		delete vec;
	}

	//値の取り出し
	// p (int*) : 得たいノイズのテクスチャ座標[0,解像度]．dim次元で指定
	// return (double) : [0,1]の値
	double at(int* p) {
		for (int i = 0; i < dim; i++)
		{
			if (p[i] < 0 || p[i] > res[i])
				throw "Texture access error!\n";
		}

		int idx = 0;
		for (int i = 0; i < dim; i++)
		{
			//次元に併せてインデックスにバイアスをかけて足す
			//e.g. (1D) x, (2D) + y * width, (3D) + z * height * width
			int index = p[i];
			for (int k = 0; k < i; k++)
			{
				index *= res[k];
			}
			idx += index;
		}
		return texture[idx];
	}

	//値の取り出し
	// p (int*) : 得たいノイズのテクスチャ座標[0,解像度]．dim次元で指定
	// return (double) : [0,1]の値
	double operator()(int* p) {
		return at(p);
	}

	//ノイズのリセット
	void reset(int seed = -1) {
		delete aexp;
		delete texture;		
		delete noise;
		init(seed);
		next();
	}

private:
	int dim;
	int* res;
	int size;
	int ts;
	double dt;
	double oct;
	double lac;
	double H;
	T* noise;
	int time;

	double* aexp;
	double* texture;

	const int INTERNAL_NOISE_RESOLUTION;

	//ノイズや内部パラメータの初期化
	void init(int seed) {
		time = 0;

		aexp = new double[(int)oct + 1];
		texture = new double[size];

		double freq = 1.0;
		for (int i = 0; i < (int)oct + 1; i++)
		{
			aexp[i] = pow(freq, -H);
			freq *= lac;
		}

		int* internal_resolution = new int[dim + 1];
		for (int i = 0; i < dim; i++)
			internal_resolution[i] = INTERNAL_NOISE_RESOLUTION;
		internal_resolution[dim] = ts;
		noise = new T(dim + 1, internal_resolution, seed);
		delete internal_resolution;
	}

	//ノイズを複数の周波数で足し合わせて座標における値を得る
	// vec (double*) : ノイズを得る座標[0,1]
	// return (double) : 実数値（オクターブ等で値の範囲は変わる）
	double get(double* vec) {
		double val = 0.0;
		int i;

		for (i = 0; i < oct; i++)
		{
			val += noise->get(vec) * aexp[i];
			for (int k = 0; k < dim + 1; k++)
				vec[k] *= lac;
		}

		double remainder = oct - (int)oct;
		if (!remainder)
			val += remainder * noise->get(vec) * aexp[i];

		return val;
	}
};