#pragma once
#include <random>

//参考：http://www.scratchapixel.com/index.php

//N次元のValue Noiseを作成するマスタークラス
class KWNoise
{
public:
	// 初期化
	// dim (uint) : ノイズの次元
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	KWNoise(unsigned int dim, int seed) :
		mt(seed),
		ureal(0.0, 1.0)
	{
		dimensions = dim;

		if (seed < 0) {
			std::random_device rdev;
			mt.seed(rdev());
		}
	}
	virtual ~KWNoise() {
	};

	// 値の取り出し
	// position (double*) : 得たいノイズの座標．dim次元で指定
	// return (double) : [0,1]の値
	virtual double get(const double* position) = 0;

	// 内部状態の初期化
	// seed (int) : 内部ノイズ生成のためのシード．-1(default)でランダム値を使用
	void reset(int seed) {
		if (seed < 0) {
			std::random_device rdev;
			mt.seed(rdev());
		}

		init();
	}

	//値の取り出し
	// position (double*) : 得たいノイズの座標．dim次元で指定
	// return (double) : [0,1]の値
	double operator()(const double* position) {
		return get(position);
	}

protected:
	//内部データの初期化
	virtual void init() = 0;
	
	//次元数
	unsigned int dimensions;

	//乱数生成器
	std::mt19937 mt;
	
	//実数均一分布
	std::uniform_real<double> ureal;

	//スムージング関数
	inline double LINEAR(double t) { return t; }
	inline double CUBIC(double t) { return t * t * (-2.0 * t + 3.0); }
	inline double QUINTIC(double t) { return t * t * t * (10.0 + t * (-15.0 + 6.0 * t)); }
	
	//アクセス不可のコンストラクタ
	KWNoise() {}
	KWNoise(const KWNoise& a) {}
};
