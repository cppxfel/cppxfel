#include <cmath>
#include <vector>
#include <iostream>
#include "parameters.h"

double gaussian_distribution(double x, double mean, double sigma, double height)
{
	double power = 0 - pow((x - mean), 2) / (2 * sigma * sigma);
	double exp = pow(M_E, power);

	exp *= height;

	return exp;
}

double sq_diff_for_num(vector<double> wavelengths,
		vector<int> refcounts, int diff_num, double mean, double sigma,
		double height)
{
	double wavelength_lower = wavelengths[diff_num];
	double wavelength_upper = wavelengths[diff_num + 1];

	double wavelength_trial = (wavelength_upper + wavelength_lower) / 2;

	double modelled_refcount = gaussian_distribution(wavelength_trial, mean,
			sigma, height);

	double true_refcount = refcounts[diff_num];

	double square_difference = pow(true_refcount - modelled_refcount, 2);

	return square_difference;
}

double find_residual(vector<double> wavelengths, vector<int> refcounts,
		int num, double mean, double sigma, double height)
{
	double residual = 0;

	for (int i = 0; i < num - 1; i++)
	{
		double addition = sq_diff_for_num(wavelengths, refcounts, i, mean, sigma,
				height);
		residual += addition;
	}

	return residual;
}

double find_mode_height(vector<double> wavelengths,
		vector<int> refcounts, int num, int *mode_num)
{
	int max_num = -1;
	double best_count = 0;

	for (int i = 0; i < num; i++)
	{
		if (refcounts[i] > best_count)
		{
			max_num = i;
			best_count = refcounts[i];
		}
	}

	*mode_num = max_num;

	return best_count;
}

void gaussian_fit(vector<double> wavelengths, vector<int> refcounts,
		int num, double *mean, double *stdev, double *score, bool print)
{
	double height = 0;
	int mode_num = 0;
	double best_mean = 0;
	double best_stdev = 0;
	double best_height = 0;
	double mean_step = 0.005;
	double stdev_step = 0.0001;
	double height_step = 5;
	int optimised_mean = 0;
	int optimised_stdev = 0;
	double best_score = 0;

	if (wavelengths.size() == 0)
		return;

	height = find_mode_height(wavelengths, refcounts, num, &mode_num);
	best_mean = wavelengths[mode_num];

	best_stdev = 0.005 * best_mean;
	int count = 0;

	while (!(optimised_mean && optimised_stdev) && count < 500)
	{
		optimised_mean = 0;
		optimised_stdev = 0;

        double mean_trials[3] = {0, 0, 0};
		double mean_scores[3] = {0, 0, 0};
		double stdev_trials[3] = {0, 0, 0};
		double stdev_scores[3] = {0, 0, 0};
		double height_trials[3] = {0, 0, 0};
		double height_scores[3] = {0, 0, 0};
		int j = 0;
		double mean_min_score = 100000;
		int mean_min_num = -1;
		double stdev_min_score = 100000;
		int stdev_min_num = -1;
		double height_min_score = 100000;
		int height_min_num = -1;

		for (double i = best_mean - mean_step; j < 3; i += mean_step)
		{
			mean_scores[j] = find_residual(wavelengths, refcounts, num, i,
					best_stdev, height);
			mean_trials[j] = i;
			j++;
		}

		for (int i = 0; i < 3; i++)
		{
			if (mean_scores[i] < mean_min_score)
			{
				mean_min_score = mean_scores[i];
				mean_min_num = i;
			}
		}

		j = 0;

		best_mean = mean_trials[mean_min_num];

		for (double i = best_stdev - stdev_step; j < 3; i += stdev_step)
		{
			stdev_scores[j] = find_residual(wavelengths, refcounts, num,
					best_mean, i, height);
			stdev_trials[j] = i;
			j++;
		}

		for (int i = 0; i < 3; i++)
			if (stdev_scores[i] < stdev_min_score)
			{
				stdev_min_score = stdev_scores[i];
				stdev_min_num = i;
			}

		best_stdev = stdev_trials[stdev_min_num];

		j = 0;

		for (double i = best_height - height_step; j < 3; i += height_step)
		{
			height_scores[j] = find_residual(wavelengths, refcounts, num,
					best_mean, best_stdev, i);
			height_trials[j] = i;
			j++;
		}

		for (int i = 0; i < 3; i++)
			if (height_scores[i] < height_min_score)
			{
				height_min_score = height_scores[i];
				height_min_num = i;
			}

		best_height = height_trials[height_min_num];

		if (mean_min_num == 1)
		{
			mean_step /= 2;
		}

		if (stdev_min_num == 1)
		{
			stdev_step /= 2;
		}

		if (height_min_num == 1)
		{
			height_step /= 2;
		}

		if (mean_step < 0.00002)
			optimised_mean = 1;

		if (stdev_step < 0.00002)
			optimised_stdev = 1;

		best_score = stdev_scores[stdev_min_num];
		count++;
	}

	if (print)
	{
		for (int i = 0; i < num - 1; i++)
		{
			if (refcounts[i] > 0)
			{
				double wavelength_lower = wavelengths[i];
				double wavelength_upper = wavelengths[i + 1];
				double wavelength_trial = (wavelength_upper + wavelength_lower)
						/ 2;

				double modelled_refcount = gaussian_distribution(
						wavelength_trial, best_mean, best_stdev, height);

				if (modelled_refcount < 0.05)
					continue;

			}
		}
	}

	double best_stderr = best_stdev / best_mean * 2;
	best_score /= height * height;

	*score = best_score;

	(*mean) = best_mean;
	(*stdev) = best_stderr;
}
