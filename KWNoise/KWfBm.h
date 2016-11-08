#pragma once
#include "KWNoise.h"

//fBm�Ńm�C�Y���ڍ׉������e�N�X�`���Ƃ��Ĉ����N���X
//template���� T : �m�C�Y�Ɏg���N���X�iKWNoise�̔h���N���X�j
//template���� N : �����m�C�Y�Ɏg���𑜓x�i�f�t�H���g128�j 
template <class T, int N = 128>
class KWfBm
{
public:
	//������
	// dimension (int) : �e�N�X�`���̎���
	// resolutions (int*) : �e�N�X�`���̉𑜓x�Ddimension�����Ŏw��
	// time_scale (int) : ���ԕ����̉𑜓x
	// delta_time (double) : next()�̃R�[�����Ƃɐi�߂鎞�ԕ�
	// octaves (double) : �I�N�^�[�u���D�d�˕�������g���̐�
	// lacunarity (double) : �󌄐�[0,1]�D�d�˕�������g���̕�
	// hurst_index (double) : �n�[�X�g�萔[0,1]�D�t���N�^�������̃p�����[�^
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
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

	//�m�C�Y�̎��Ԃ�i�߂�
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

	//�l�̎��o��
	// p (int*) : �������m�C�Y�̃e�N�X�`�����W[0,�𑜓x]�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	double at(int* p) {
		for (int i = 0; i < dim; i++)
		{
			if (p[i] < 0 || p[i] > res[i])
				throw "Texture access error!\n";
		}

		int idx = 0;
		for (int i = 0; i < dim; i++)
		{
			//�����ɕ����ăC���f�b�N�X�Ƀo�C�A�X�������đ���
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

	//�l�̎��o��
	// p (int*) : �������m�C�Y�̃e�N�X�`�����W[0,�𑜓x]�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	double operator()(int* p) {
		return at(p);
	}

	//�m�C�Y�̃��Z�b�g
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

	//�m�C�Y������p�����[�^�̏�����
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

	//�m�C�Y�𕡐��̎��g���ő������킹�č��W�ɂ�����l�𓾂�
	// vec (double*) : �m�C�Y�𓾂���W[0,1]
	// return (double) : �����l�i�I�N�^�[�u���Œl�͈͕̔͂ς��j
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