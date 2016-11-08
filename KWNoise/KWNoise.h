#pragma once
#include <random>

//�Q�l�Fhttp://www.scratchapixel.com/index.php

//N������Value Noise���쐬����}�X�^�[�N���X
class KWNoise
{
public:
	// ������
	// dim (uint) : �m�C�Y�̎���
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
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

	// �l�̎��o��
	// position (double*) : �������m�C�Y�̍��W�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	virtual double get(const double* position) = 0;

	// ������Ԃ̏�����
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
	void reset(int seed) {
		if (seed < 0) {
			std::random_device rdev;
			mt.seed(rdev());
		}

		init();
	}

	//�l�̎��o��
	// position (double*) : �������m�C�Y�̍��W�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	double operator()(const double* position) {
		return get(position);
	}

protected:
	//�����f�[�^�̏�����
	virtual void init() = 0;
	
	//������
	unsigned int dimensions;

	//����������
	std::mt19937 mt;
	
	//�����ψꕪ�z
	std::uniform_real<double> ureal;

	//�X���[�W���O�֐�
	inline double LINEAR(double t) { return t; }
	inline double CUBIC(double t) { return t * t * (-2.0 * t + 3.0); }
	inline double QUINTIC(double t) { return t * t * t * (10.0 + t * (-15.0 + 6.0 * t)); }
	
	//�A�N�Z�X�s�̃R���X�g���N�^
	KWNoise() {}
	KWNoise(const KWNoise& a) {}
};
