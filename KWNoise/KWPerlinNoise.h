#pragma once
#include "KWNoise.h"
#include <algorithm>

//�I���W�i���i���܂�Q�Ƃ��ĂȂ��j
//http://mrl.nyu.edu/~perlin/noise/

//N������Perlin Noise���쐬����N���X
//N�����̓��������_��N�����x�N�g���������ߗe�ʂ��傫��
//�����m�C�Y�̉𑜓x�̓��[�U��N�����Ŏw�肷��
//�l�̎��o����[0,1]�͈͂̍��W�w��ōs��
//[0,1]�𒴂���ۂ͏��������g���̂ŌJ��Ԃ����m�C�Y�����������
//�������傫������ƈ����Ȃ��D32bit�����̂��ߓ����L���̈�𖳎����Ă�32���܂ł�
class KWPerlinNoise : public KWNoise
{
public:
	// ������
	// dim (uint) : �m�C�Y�̎���
	// resolution (int*) : �����m�C�Y�̉𑜓x�Ddim�����Ŏw��
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
	KWPerlinNoise(unsigned int dim, const int* resolutions, int seed = -1) :
		KWNoise(dim, seed)
	{
		//�𑜓x�̐ݒ�
		res = new int[dimensions];
		size = 1;
		for (unsigned int i = 0; i < dimensions; i++)
		{
			size *= (res[i] = resolutions[i]);
		}
		data = new double[size * dimensions];

		//�����m�C�Y�̏�����
		init();

		//�v�Z�̈�̊m��
		index = new int[2 * dimensions];
		remainders = new double[dimensions];
		vec_p = new double[2 * dimensions];
		v0id = new int[dimensions];
		v1id = new int[dimensions];
	}

	~KWPerlinNoise() {
		delete index;
		delete remainders;
		delete vec_p;
		delete v0id;
		delete v1id;
		delete res;
		delete data;
	}

	// �l�̎��o��
	// position (double*) : �������m�C�Y�̍��W�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	double get(const double* position) {

		//���W�̏���
		for (unsigned int i = 0; i < dimensions; i++)
		{
			//[0,1]�̒l��[0,�𑜓x-1]�ɒ����D�������̂ݎg��
			double int_part;
			double p = std::modf(position[i], &int_part) * (res[i] - 1);

			// (1D)
			//   <-- t -->
			//  idx0      p    idx1
			// -|---------+----|-
			//
			//�����W���ł��߂��O��2�̐������W�ɒ����]�蕪�����߂�
			int idx0 = (int)std::floor(p); //�O���W
			int idx1 = (idx0 < res[i] - 1) ? idx0 + 1 : idx0; //����W
			double t = p - (int)std::floor(p);
			
			//�����ɕ����ăC���f�b�N�X�Ƀo�C�A�X��������ɉ��Z���邾���ň�����悤�ɂ���
			//e.g. (1D) x, (2D) y * width, (3D) z * height * width
			for (unsigned int k = 0; k < i; k++)
			{
				idx0 *= res[i];
				idx1 *= res[i];
			}

			//�C���f�b�N�X�Ɨ]��l���i�[
			//index��1D�̃C���f�b�N�X0,1�C2D�̃C���f�b�N�X0,1�C...�Ɗi�[
			index[i * 2 + 0] = idx0;
			index[i * 2 + 1] = idx1;
			remainders[i] = CUBIC(t); //�]��l��5���֐��̃X���[�W���O�W���ɂ���

			// (1D)
			//    P0      p  P1
			// -|-------->+<---|-
			//�i�q�_������W�_�ւ̃x�N�g���iP�x�N�g���j���L�^����
			vec_p[i * 2 + 0] = t;
			vec_p[i * 2 + 1] = t - 1.0;
		}

		//1D�ɂ������ԏ���
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
		// �܂���1�����ł̕�ԁC�܂�o�ʒu�ł̒l�����߂ăo�b�t�@�Ɋi�[
		//
		//���W�̑g�ݍ��킹��N�����ł�2^N����
		//�����̏�����1������2�g�݂̒l�͕�Ԃ���̂�2^(N-1)�̒l�̕ۑ��p�Ƀo�b�t�@���m��
		int n = (int)pow(2, dimensions - 1);
		double* buff = new double[n];
		//1D�ł̍��W�C���f�b�N�X
		int idx0 = index[0];
		int idx1 = index[1];
		//1D�ł�P�x�N�g���C���f�b�N�X
		v0id[0] = 0;
		v1id[0] = 1;
		//2^N�̍��W����������D1���[�v��1������2�g�̍��W���g���̂Ń��[�v��2^(N-1)
		for (int b = 0; b < n; b++)
		{
			//2�����ȏ�ň������W�̑g�ݍ��킹���r�b�g�ōl����
			//e.g. (3D) b=2_(10) = 10_(2) => y1, z0

			//�������W�̃C���f�b�N�X��P�x�N�g���̃C���f�b�N�X�����߂�
			int idx = 0;
			for (unsigned int k = 0; k < dimensions - 1; k++)
			{
				//�E����k�Ԗڂ̃r�b�g�������Ă���Ό���W���g��
				if (b & (1 << k)) {
					idx += index[2 * (k + 1) + 1];
					v0id[k + 1] = v1id[k + 1] = 2 * (k + 1) + 1;
				}
				else {
					idx += index[2 * (k + 1) + 0];
					v0id[k + 1] = v1id[k + 1] = 2 * (k + 1) + 0;
				}
			}

			//���W�ł̃x�N�g���iC�x�N�g���j��P�x�N�g���̓��ς�����Ă��̈ʒu�̒l�Ƃ���
			// (2D)
			//            x1,y0
			// -x---------o----x-
			//  |              |
			//  |         *    |
			//  |         p    |
			//  |              |
			// �����Ő�������x��̃x�N�g����P�x�N�g���ix����*�ւ̃x�N�g���j��x�ł̒l�Ƃ��č̗p
			// 1D�����ŕ�Ԃ���o��̒l���o�b�t�@�ɕۑ�����
			//
			double dot_product0 = 0.0;
			double dot_product1 = 0.0;
			for (unsigned int k = 0; k < dimensions; k++)
			{
				dot_product0 += data[(idx + idx0) * dimensions + k] * (vec_p[v0id[k]]);
				dot_product1 += data[(idx + idx1) * dimensions + k] * (vec_p[v1id[k]]);
			}
			buff[b] = dot_product0 +remainders[0] * (dot_product1 - dot_product0);
		}

		//2�����ȍ~�̕�ԏ���
		//���łɃo�b�t�@��1�����ŕ�Ԃ����l�����Ɋi�[����Ă���
		//�O����2�����Ԃ��Ă����Ηǂ�
		for (unsigned int i = 1; i < dimensions; i++)
		{
			//�O����2����i�����ڂ̗]��l�ŕ�Ԃ��Ă��̒l���o�b�t�@�̐擪�ɋl�߂Ă���
			int cnt = 0;//�l�߂Ă����o�b�t�@�̐擪�C���f�b�N�X
			for (int k = 0; k < n; k += 2) {
				double a = buff[k + 0];
				buff[k + 0] = 0.0;
				double b = buff[k + 1];
				buff[k + 1] = 0.0;
				buff[cnt++] = a + remainders[i] * (b - a);
			}
			n /= 2;
		}

		//�o�b�t�@�̐擪���ŏI�I�ɕ�Ԃ��ꂽ�l
		//�o�b�t�@�̉��������̂ŕϐ��ɑҔ�
		double ret = buff[0];
		delete buff;

		//���ς̒l��[-1,1]�Ȃ̂�Value Noise�Ǝd�l�𕹂��邽��[0,1]�ɒ���
		return (ret + 1.0) / 2.0;
	}

	// ������Ԃ̏�����
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
	void reset(int seed = -1) {
		KWNoise::reset(seed);
	}

protected:
	void init() {
		for (int i = 0; i < size; i++)
		{
			//�����m�C�Y��[0,1]�̗����ŏ���������

			//N�����̃����_���ȃx�N�g���Ƃ��Đ��K�����ꂷ�邽�ߒ����𓾂�
			double length = 0.0;
			for (unsigned int k = 0; k < dimensions; k++)
			{
				double v = data[i * dimensions + k] = 2.0 * ureal(mt) - 1.0;
				length += v * v;
			}
			//���K������D����̃x�N�g���v�f��data���ɘA�������ʒu�ɕۑ�����
			length = sqrt(length);
			for (unsigned int k = 0; k < dimensions; k++)
			{
				data[i * dimensions + k] /= length;
			}
		}
	}

private:
	int* res;
	int size;
	double* data;
	int* index;
	double* remainders;
	double* vec_p;
	int* v0id;
	int* v1id;
};


//N������Perlin Noise���쐬����N���X
//�J�ڃe�[�u����p���ē����m�C�Y���y�ʉ����Ă���
//�����m�C�Y�̉𑜓x�̓��[�U���w�肷�邪1������2�̙p�搔�Ɍ��肳���
//�l�̎��o���͔C�ӂ̐��ōs���邪
//�e�N�X�`��������ۂɂ͏����_�ȉ��̉𑜓x�����������W�������K�͂ōL������K�v������
//i.e [0,1]���ƃm�C�Y�����ꂸ�C�������݂��ƃz���C�g�m�C�Y�ɂȂ�
//�������傫������ƈ����Ȃ��D32bit�����̂��ߓ����m�C�Y�T�C�Y�𖳎����Ă�32���܂ł�
class KWPerlinNoiseLight : public KWNoise
{
public:
	// ������
	// dim (uint) : �m�C�Y�̎���
	// resolution (int*) : �����m�C�Y�̉𑜓x�Ddim�����Ŏw��
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
	KWPerlinNoiseLight(unsigned int dim, const int* resolution, int seed = -1) :
		KWNoise(dim, seed)
	{
		//�𑜓x�̐ݒ�
		res = resolution[0];
		data = new double[res * dimensions];
		permutation_table = new int[2 * res];

		//�����m�C�Y�̏�����
		init();

		//�v�Z�̈�̊m��
		index = new int[2 * dimensions];
		remainders = new double[dimensions];
		vec_p = new double[2 * dimensions];
		v0id = new int[dimensions];
		v1id = new int[dimensions];
		idx0 = new int[dimensions];
		idx1 = new int[dimensions];
	}

	~KWPerlinNoiseLight() {
		delete idx0;
		delete idx1;
		delete index;
		delete remainders;
		delete vec_p;
		delete v0id;
		delete v1id;
		delete permutation_table;
		delete data;
	}

	// �l�̎��o��
	// position (double*) : �������m�C�Y�̍��W�Ddim�����Ŏw��
	// return (double) : [0,1]�̒l
	double get(const double* position) {

		//���W�̏���
		for (unsigned int i = 0; i < dimensions; i++)
		{
			double p = position[i];

			// (1D)
			//   <-- t -->
			//  idx0      p    idx1
			// -|---------+----|-
			//
			//[0,res-1]�͈͂ɃN�����v��
			//�����W���ł��߂��O��2�̐������W�ɒ����]�蕪�����߂�
			int idx0 = (int)std::floor(p) & (res - 1); //�O���W
			int idx1 = (idx0 + 1) & (res - 1); //����W
			double t = p - (int)std::floor(p);
			
			//�C���f�b�N�X�Ɨ]��l���i�[
			//index��1D�̃C���f�b�N�X0,1�C2D�̃C���f�b�N�X0,1�C...�Ɗi�[
			index[i * 2 + 0] = idx0;
			index[i * 2 + 1] = idx1;
			remainders[i] = CUBIC(t); //�]��l��5���֐��̃X���[�W���O�W���ɂ���
			
			// (1D)
			//    P0      p  P1
			// -|-------->+<---|-
			//�i�q�_������W�_�ւ̃x�N�g���iP�x�N�g���j���L�^����
			vec_p[i * 2 + 0] = t; //P0
			vec_p[i * 2 + 1] = t - 1.0; //P1
		}

		//1D�ɂ������ԏ���
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
		// �܂���1�����ł̕�ԁC�܂�o�ʒu�ł̒l�����߂ăo�b�t�@�Ɋi�[
		//
		//���W�̑g�ݍ��킹��N�����ł�2^N����
		//�����̏�����1������2�g�݂̒l�͕�Ԃ���̂�2^(N-1)�̒l�̕ۑ��p�Ƀo�b�t�@���m��
		int n = (int)pow(2, dimensions - 1);
		double* buff = new double[n];
		//1D�ł̍��W
		idx0[0] = index[0];
		idx1[0] = index[1];
		//1D�ł�P�x�N�g���C���f�b�N�X
		v0id[0] = 0;
		v1id[0] = 1;
		//2^N�̍��W����������D1���[�v��1������2�g�̍��W���g���̂Ń��[�v��2^(N-1)
		for (int b = 0; b < n; b++)
		{
			//2�����ȏ�ň������W�̑g�ݍ��킹���r�b�g�ōl����
			//e.g. (3D) b=2_(10) = 10_(2) => y1, z0

			//�������W�̃C���f�b�N�X��P�x�N�g���̃C���f�b�N�X�����߂�
			for (unsigned int k = 0; k < dimensions - 1; k++)
			{
				//�E����k�Ԗڂ̃r�b�g�������Ă���Ό���W���g��
				if (b & (1 << k)) {
					idx0[k + 1] = idx1[k + 1] = index[2 * (k + 1) + 1];
					v0id[k + 1] = v1id[k + 1] = 2 * (k + 1) + 1;
				}
				else {
					idx0[k + 1] = idx1[k + 1] = index[2 * (k + 1) + 0];
					v0id[k + 1] = v1id[k + 1] = 2 * (k + 1) + 0;
				}
			}

			//�����̍��W�l��p�����n�b�V���֐��ɒʂ��C���f�b�N�X�𓾂�
			//�n�b�V���֐��͍��W�ɂ��Ĉ�ӂȒl��^���邽�߂��̍��W�ł̃m�C�Y�l�����߂�C���f�b�N�X�Ƃ��ē���
			int i0 = hash(idx0);
			int i1 = hash(idx1);

			//���W�ł̃x�N�g���iC�x�N�g���j��P�x�N�g���̓��ς�����Ă��̈ʒu�̒l�Ƃ���
			// (2D)
			//            x1,y0
			// -x---------o----x-
			//  |              |
			//  |         *    |
			//  |         p    |
			//  |              |
			// �����Ő�������x��̃x�N�g����P�x�N�g���ix����*�ւ̃x�N�g���j��x�ł̒l�Ƃ��č̗p
			// 1D�����ŕ�Ԃ���o��̒l���o�b�t�@�ɕۑ�����
			//
			double dot_product0 = 0.0;
			double dot_product1 = 0.0;
			for (unsigned int k = 0; k < dimensions; k++)
			{
				dot_product0 += data[i0 * dimensions + k] * (vec_p[v0id[k]]);
				dot_product1 += data[i1 * dimensions + k] * (vec_p[v1id[k]]);
			}
			buff[b] = dot_product0 + remainders[0] * (dot_product1 - dot_product0);
		}

		//2�����ȍ~�̕�ԏ���
		//���łɃo�b�t�@��1�����ŕ�Ԃ����l�����Ɋi�[����Ă���
		//�O����2�����Ԃ��Ă����Ηǂ�
		for (unsigned int i = 1; i < dimensions; i++)
		{
			//�O����2����i�����ڂ̗]��l�ŕ�Ԃ��Ă��̒l���o�b�t�@�̐擪�ɋl�߂Ă���
			int cnt = 0;//�l�߂Ă����o�b�t�@�̐擪�C���f�b�N�X
			for (int k = 0; k < n; k += 2) {
				double a = buff[k + 0];
				buff[k + 0] = 0.0;
				double b = buff[k + 1];
				buff[k + 1] = 0.0;
				buff[cnt++] = a + remainders[i] * (b - a);
			}
			n /= 2;
		}

		//�o�b�t�@�̐擪���ŏI�I�ɕ�Ԃ��ꂽ�l
		//�o�b�t�@�̉��������̂ŕϐ��ɑҔ�
		double ret = buff[0];
		delete buff;

		//���ς̒l��[-1,1]�Ȃ̂�Value Noise�Ǝd�l�𕹂��邽��[0,1]�ɒ���
		return (ret + 1.0) / 2.0;
	}

	// ������Ԃ̏�����
	// seed (int) : �����m�C�Y�����̂��߂̃V�[�h�D-1(default)�Ń����_���l���g�p
	void reset(int seed = -1) {
		KWNoise::reset(seed);
	}

protected:
	void init() {
		//�����m�C�Y��[0,1]�̗����ŏ���������
		//�����đJ�ڃe�[�u�������
		//�J�ڃe�[�u����res-1�̑S�Ă̒l������
		for (int i = 0; i < res; i++)
		{
			//�����m�C�Y��[-1,1]�̗����ŏ���������
			//N�����̃����_���ȃx�N�g���Ƃ��Đ��K�����ꂷ�邽�ߒ����𓾂�
			double length = 0.0;
			for (unsigned int k = 0; k < dimensions; k++)
			{
				double v = data[i * dimensions + k] = 2.0 * ureal(mt) - 1.0;
				length += v * v;
			}
			//���K������D����̃x�N�g���v�f��data���ɘA�������ʒu�ɕۑ�����
			length = sqrt(length);
			for (unsigned int k = 0; k < dimensions; k++)
			{
				data[i * dimensions + k] /= length;
			}
			permutation_table[i] = i;
		}

		//�J�ڃe�[�u���̒l���V���b�t������
		//�e�[�u����res�Ԗڈȍ~��0����̌J��Ԃ�
		std::uniform_int_distribution<> uinteger(0, res - 1);
		for (int i = 0; i < res; i++)
		{
			//[0,res-1]�̗���idx�Ԗڂ�i�Ԗڂ��X���b�v����
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
	double* vec_p;
	int* v0id;
	int* v1id;

	//N�����̐������W����J�ڃe�[�u�����̐����l�����߂�
	int hash(const int* p) {
		int a = 0;
		for (unsigned int i = 0; i < dimensions; i++)
			a = permutation_table[a + p[i]];
		return a;
	}
};