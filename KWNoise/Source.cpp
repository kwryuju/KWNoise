#include "KWValueNoise.h"
#include "KWPerlinNoise.h"
#include "KWfBm.h"
#include <iostream>
#include "opencv2\core.hpp"
#include "opencv2\highgui.hpp"

int main() {
	
	int SEED = 6;

	int W = 512;
	int H = 512;
	int D = 128;
	
	int RES[] = { 64, 64, 64 };
	int res[] = {W, H};

	//Value Noiseのインスタンシング
	cv::Mat val_img(H, W, CV_8UC3, cv::Scalar(0.0, 0.0, 0.0));
	KWNoise* vnoise = new KWValueNoise(3, RES, SEED);
	//KWNoise* vnoise = new KWValueNoiseLight(3, RES, SEED);
	cv::namedWindow("value noise (raw)", cv::WINDOW_KEEPRATIO);

	//Perlin Noiseのインスタンシング
	cv::Mat per_img(H, W, CV_8UC3, cv::Scalar(0.0, 0.0, 0.0));
	KWNoise* pnoise = new KWPerlinNoise(3, RES, SEED);
	//KWNoise* pnoise = new KWPerlinNoiseLight(3, RES, SEED);
	cv::namedWindow("Perlin noise (raw)", cv::WINDOW_KEEPRATIO);
	
	//Value NoiseでのfBmのインスタンシング
	cv::Mat vfBm_img(H, W, CV_8UC3, cv::Scalar(0.0, 0.0, 0.0));
	KWfBm<KWValueNoise> vfBm(2, res, D, 1.0, 4.0, 0.5, 1.0, SEED);
	cv::namedWindow("value noise (fBm)", cv::WINDOW_KEEPRATIO);

	//Perlin NoiseでのfBmのインスタンシング
	cv::Mat pfBm_img(H, W, CV_8UC3, cv::Scalar(0.0, 0.0, 0.0));
	KWfBm<KWPerlinNoise> pfBm(2, res, D, 1.0, 4.0, 0.5, 1.0, SEED);
	cv::namedWindow("Perlin noise (fBm)", cv::WINDOW_KEEPRATIO);

	//ノイズ生成ループ（any keyで時間を進める）
	int z = 0;
	while (true)
	{
		for (int y = 0; y < H; y++)
		{
			for (int x = 0; x < W; x++)
			{
				double scale_x = 1.0f / (double)(W - 1);
				double scale_y = 1.0f / (double)(H - 1);
				double scale_z = 1.0f / (double)(D - 1);
				double p0[] = { x * scale_x, y * scale_y, z * scale_z };
				int p1[] = { x, y };

				//Value Noise
				double val = vnoise->get(p0);
				val_img.at<cv::Vec3b>(y, x)[0] =
					val_img.at<cv::Vec3b>(y, x)[1] =
					val_img.at<cv::Vec3b>(y, x)[2] =
					(uchar)(val * 255);

				//Perlin Noise
				val = pnoise->get(p0);
				per_img.at<cv::Vec3b>(y, x)[0] =
					per_img.at<cv::Vec3b>(y, x)[1] =
					per_img.at<cv::Vec3b>(y, x)[2] =
					(uchar)(val * 255);

				//fBm (Value Noise)
				val = vfBm(p1);
				vfBm_img.at<cv::Vec3b>(y, x)[0] =
					vfBm_img.at<cv::Vec3b>(y, x)[1] =
					vfBm_img.at<cv::Vec3b>(y, x)[2] =
					(uchar)(val * 255);

				//fBm (Perlin Noise)
				val = pfBm(p1);
				pfBm_img.at<cv::Vec3b>(y, x)[0] =
					pfBm_img.at<cv::Vec3b>(y, x)[1] =
					pfBm_img.at<cv::Vec3b>(y, x)[2] =
					(uchar)(val * 255);
			}
		}

		//時間を進める
		z = (z < D) ? z + 1 : 0;
		vfBm.next();
		pfBm.next();

		//表示
		cv::imshow("value noise (raw)", val_img);
		cv::imshow("Perlin noise (raw)", per_img);
		cv::imshow("value noise (fBm)", vfBm_img);
		cv::imshow("Perlin noise (fBm)", pfBm_img);

		//入力待ち
		auto c = cv::waitKey(0);
		if (c == 13) break; //Enterで終了
	}


	//解放処理
	cv::destroyAllWindows();
	delete vnoise;
	delete pnoise;

	return 0;
}