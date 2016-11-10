#pragma once
#include "KWNoise.h"

//N次元のValue Noiseを作成するクラス
//N次元の内部ノイズ(WhiteNoise)を持つため容量が大きい
//内部ノイズの解像度はユーザがN次元で指定する
//値の取り出しは[0,1]範囲の座標指定で行う
//[0,1]を超える際は小数部を使うので繰り返しや非ノイズが生成される
//次元が大きすぎると扱えない．32bit実装のため内部ノイズサイズを無視しても32次までか
class KWValueNoise : public KWNoise
{
public:
	// 初期化
	// dim (uint) : ノイズの次元
	// resolution (int*) : 内部ノイズの解像度．dim次元で指定
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	KWValueNoise(unsigned int dim, const int* resolutions, int seed = -1) :
		KWNoise(dim, seed)
	{
		//解像度の設定
		res = new int[dimensions];
		size = 1;
		for (unsigned int i = 0; i < dimensions; i++)
			size *= (res[i] = resolutions[i]);
		data = new double[size];

		//内部ノイズの初期化
		init();

		//計算領域の確保
		index = new int[2 * dimensions];
		remainders = new double[dimensions];
	}

	~KWValueNoise() {
		delete[] index;
		delete[] remainders;
		delete[] res;
		delete[] data;
	}

	// 値の取り出し
	// position (double*) : 得たいノイズの座標．dim次元で指定
	// return (double) : [0,1]の値
	double get(const double* position) {
		
		//座標の処理
		for (unsigned int i = 0; i < dimensions; i++)
		{
			//[0,1]の値を[0,解像度-1]に直す．小数部のみ使う
			double int_part;
			double p = std::modf(position[i], &int_part) * (res[i] - 1);

			// (1D)
			//   <-- t -->
			//  idx0      p    idx1
			// -|---------+----|-
			//
			//実座標を最も近い前後2つの整数座標に直し余剰分を求める
			int idx0 = (int)std::floor(p); //前座標
			int idx1 = (idx0 < res[i] - 1) ? idx0 + 1 : idx0; //後座標
			double t = p - (int)std::floor(p);

			//次元に併せてインデックスにバイアスをかけ後に加算するだけで扱えるようにする
			//e.g. (1D) x, (2D) y * width, (3D) z * height * width
			for (unsigned int k = 0; k < i; k++)
			{
				idx0 *= res[i];
				idx1 *= res[i];
			}

			//インデックスと余剰値を格納
			//indexは1Dのインデックス0,1，2Dのインデックス0,1，...と格納
			index[i * 2 + 0] = idx0;
			index[i * 2 + 1] = idx1;
			remainders[i] = CUBIC(t); //余剰値を5次関数のスムージング係数にする
		}

		//1Dにおける補間処理
		//
		// (2D)
		//  x0,y0          x1,y0
		// -|---------o----|-
		//  |              |
		//  |              |
		//  |         p    |
		//  +         *    +
		//  |              |
		// -|---------o----|-
		//  x0,y1          x1,y1
		// まずは1次元での補間，つまりo位置での値を求めてバッファに格納
		//
		//座標の組み合わせはN次元では2^N個ある
		//ここの処理で1次元の2個組みの値は補間するので2^(N-1)個の値の保存用にバッファを確保
		int n = (int)pow(2, dimensions - 1);
		double* buff = new double[n];
		//1Dでの座標インデックス
		int idx0 = index[0];
		int idx1 = index[1];
		//2^N個の座標を処理する．1ループで1次元で2個組の座標を使うのでループは2^(N-1)
		for (int b = 0; b < n; b++)
		{
			//2次元以上で扱う座標の組み合わせをビットで考える
			//e.g. (3D) b=2_(10) = 10_(2) => y1, z0

			//扱う座標のインデックスを求める
			int idx = 0;
			for (unsigned int k = 0; k < dimensions - 1; k++)
			{
				//右からk番目のビットが立っていれば後座標を使う
				if (b & (1 << k))
					idx += index[2 * (k + 1) + 1];
				else
					idx += index[2 * (k + 1) + 0];
			}

			//1次元方向で補間する
			//e.g. (3D) b=2 => (x0, y1, z0)と(x1, y1, z0)を1次元の余剰値で補間
			buff[b] = data[idx + idx0] + remainders[0] * (data[idx + idx1] - data[idx + idx0]);
		}
		
		//2次元以降の補間処理
		//すでにバッファに1次元で補間した値が順に格納されている
		//前から2つずつを補間していけば良い
		for (unsigned int i = 1; i < dimensions; i++)
		{
			//前から2つずつi次元目の余剰値で補間してその値をバッファの先頭に詰めていく
			int cnt = 0;//詰めていくバッファの先頭インデックス
			for (int k = 0; k < n; k += 2) {
				double a = buff[k + 0];
				buff[k + 0] = 0.0;
				double b = buff[k + 1];
				buff[k + 1] = 0.0;
				buff[cnt++] = a + remainders[i] * (b - a);
			}
			n /= 2;
		}

		//バッファの先頭が最終的に補間された値
		//バッファの解放があるので変数に待避
		double ret = buff[0];
		delete[] buff;

		return ret;
	}

	// 内部状態の初期化
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	void reset(int seed = -1) {
		KWNoise::reset(seed);
	}

protected:
	void init() {
		//内部ノイズを[0,1]の乱数で初期化する
		for (int i = 0; i < size; i++)
		{
			data[i] = ureal(mt);
		}
	}

private:
	int* res;
	int size;
	double* data;
	int* index;
	double* remainders;
};


//N次元のValue Noiseを作成するクラス
//遷移テーブルを用いて内部ノイズを軽量化している
//内部ノイズの解像度はユーザが指定するが1次元で2の冪乗数に限定される
//値の取り出しは任意の数で行えるが
//テクスチャ化する際には小数点以下の解像度を持ちレンジが整数規模で広くある必要がある
//i.e [0,1]だとノイズ化されず，整数刻みだとホワイトノイズになる
//次元が大きすぎると扱えない．32bit実装のため内部ノイズサイズを無視しても32次までか
class KWValueNoiseLight : public KWNoise
{
public:
	// 初期化
	// dim (uint) : ノイズの次元
	// resolution (int*) : 内部ノイズの解像度．resolution[0]のみ使用
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	KWValueNoiseLight(unsigned int dim, const int* resolution, int seed = -1) :
		KWNoise(dim, seed)
	{
		//解像度の設定
		res = resolution[0];
		data = new double[res];
		permutation_table = new int[2 * res];

		//内部ノイズの初期化
		init();

		//計算領域の確保
		index = new int[2 * dimensions];
		remainders = new double[dimensions];
		idx0 = new int[dimensions];
		idx1 = new int[dimensions];
	}

	~KWValueNoiseLight() {
		delete[] idx0;
		delete[] idx1;
		delete[] index;
		delete[] remainders;
		delete[] permutation_table;
		delete[] data;
	}

	// 値の取り出し
	// position (double*) : 得たいノイズの座標．dim次元で指定
	// return (double) : [0,1]の値
	double get(const double* position) {

		//座標の処理
		for (unsigned int i = 0; i < dimensions; i++)
		{
			double p = position[i];
			
			// (1D)
			//   <-- t -->
			//  idx0      p    idx1
			// -|---------+----|-
			//
			//[0,res-1]範囲にクランプし
			//実座標を最も近い前後2つの整数座標に直し余剰分を求める
			int idx0 = (int)(std::floor(p)) & (res - 1); //前座標
			int idx1 = (idx0 + 1) & (res - 1); //後座標
			double t = p - (int)(std::floor(p));//idx0はクランプされているので使わない
			
			//インデックスと余剰値を格納
			//indexは1Dのインデックス0,1，2Dのインデックス0,1，...と格納
			index[2 * i + 0] = idx0;
			index[2 * i + 1] = idx1;
			remainders[i] = CUBIC(t); //余剰値を5次関数のスムージング係数にする
		}

		//1Dにおける補間処理
		//
		// (2D)
		//  x0,y0          x1,y0
		// -|---------o----|-
		//  |              |
		//  |              |
		//  |         p    |
		//  +         *    +
		//  |              |
		// -|---------o----|-
		//  x0,y1          x1,y1
		// まずは1次元での補間，つまりo位置での値を求めてバッファに格納
		//
		//座標の組み合わせはN次元では2^N個ある
		//ここの処理で1次元の2個組みの値は補間するので2^(N-1)個の値の保存用にバッファを確保
		int n = (int)pow(2, dimensions - 1);
		double* buff = new double[n];
		//1Dでの座標
		idx0[0] = index[0];
		idx1[0] = index[1];
		//2^N個の座標を処理する．1ループで1次元で2個組の座標を使うのでループは2^(N-1)
		for (int b = 0; b < n; b++)
		{
			//2次元以上で扱う座標の組み合わせをビットで考える
			//e.g. (3D) b=2_(10) = 10_(2) => y1, z0

			//扱う座標を求める
			for (unsigned int k = 0; k < dimensions - 1; k++)
			{
				//右からk番目のビットが立っていれば後座標を使う
				if (b & (1 << k))
					idx0[k + 1] = idx1[k + 1] = index[2 * (k + 1) + 1];
				else
					idx0[k + 1] = idx1[k + 1] = index[2 * (k + 1) + 0];
			}

			//次元個の座標値を用いたハッシュ関数に通しインデックスを得る
			//ハッシュ関数は座標について一意な値を与えるためその座標でのノイズ値を決めるインデックスとして働く
			int i0 = hash(idx0);
			int i1 = hash(idx1);

			//1次元方向で補間する
			//e.g. (3D) b=2 => (x0, y1, z0)と(x1, y1, z0)を1次元の余剰値で補間
			buff[b] = data[i0] + remainders[0] * (data[i1] - data[i0]);
		}

		//2次元以降の補間処理
		//すでにバッファに1次元で補間した値が順に格納されている
		//前から2つずつを補間していけば良い
		for (unsigned int i = 1; i < dimensions; i++)
		{
			//前から2つずつi次元目の余剰値で補間してその値をバッファの先頭に詰めていく
			int cnt = 0;//詰めていくバッファの先頭インデックス
			for (int k = 0; k < n; k += 2) {
				double a = buff[k + 0];
				buff[k + 0] = 0.0;
				double b = buff[k + 1];
				buff[k + 1] = 0.0;
				buff[cnt++] = a + remainders[i] * (b - a);
			}
			n /= 2;
		}

		//バッファの先頭が最終的に補間された値
		//バッファの解放があるので変数に待避
		double ret = buff[0];
		delete[] buff;

		return ret;
	}

	// 内部状態の初期化
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	void reset(int seed = -1) {
		KWNoise::reset(seed);
	}

protected:

	void init() {
		//内部ノイズを[0,1]の乱数で初期化する
		//加えて遷移テーブルも作る
		//遷移テーブルはres-1の全ての値を持つ
		for (int i = 0; i < res; i++)
		{
			data[i] = ureal(mt);
			permutation_table[i] = i;
		}

		//遷移テーブルの値をシャッフルする
		//テーブルのres番目以降は0からの繰り返し
		std::uniform_int_distribution<> uinteger(0, res - 1);
		for (int i = 0; i < res; i++)
		{
			//[0,res-1]の乱数idx番目とi番目をスワップする
			int idx = uinteger(mt);
			int t = permutation_table[i];
			permutation_table[res + i] = 
			permutation_table[i] = permutation_table[idx];
			permutation_table[idx] = t;
		}
	}

private:
	int res;
	double* data;
	int* permutation_table;
	int* index;
	int* idx0;
	int* idx1;
	double* remainders;

	//N次元の整数座標から遷移テーブル内の整数値を求める
	int hash(const int* p) {
		int a = 0;
		for (unsigned int i = 0; i < dimensions; i++)
			a = permutation_table[a + p[i]];
		return a;
	}
};