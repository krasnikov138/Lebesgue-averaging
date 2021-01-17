#include <vector>

#include "endf.h"
#include "spectrum.h"

#ifndef DISTRIBUTION
#define DISTRIBUTION

struct Distribution2D {
	size_t dim;
    LinearGrid primary_energies;
    std::vector<LinearGrid> secondary_energies;
    std::vector<std::vector<std::vector<double>>> legandre_coefs;

    Distribution2D(): dim(0) {};
    Distribution2D(const ContinuumDistribution& distr, int dim = -1);

    size_t max_dim(const ContinuumDistribution& distr) const;
    std::vector<double> unify_legandre_polinom(const std::vector<double>& polinom) const;

    std::vector<double> interpolate(int i, double point) const;
    std::vector<double> interpolate(double primary_point, double point1, double point2) const;
};

size_t Distribution2D::max_dim(const ContinuumDistribution& distr) const {
	size_t maxd = 0;
	const auto& coefs = distr.coefs;
	for (size_t i = 0; i < coefs.size(); ++i)
		for (size_t j = 0; j < coefs[i].size(); ++j)
			maxd = std::max(maxd, coefs[i][j].size());
	return maxd;
}

std::vector<double> Distribution2D::unify_legandre_polinom(const std::vector<double>& polinom) const {
	std::vector<double> result(dim);
	size_t copy_dim = std::min(polinom.size(), dim);
	std::copy(polinom.cbegin(), polinom.cbegin() + copy_dim, result.begin());
	return result;
}

Distribution2D::Distribution2D(const ContinuumDistribution& distr, int dim) {
	primary_energies = distr.primary_energies;
	secondary_energies.reserve(distr.secondary_energies.size());
	for (const auto& energies: distr.secondary_energies)
		secondary_energies.push_back(LinearGrid(energies));

	this->dim = (dim < 0) ? max_dim(distr) : dim;

	const auto& coefs = distr.coefs;
	legandre_coefs.reserve(coefs.size());
	for (size_t i = 0; i < coefs.size(); ++i) {
		std::vector<std::vector<double>> secondary_dist;
		for (size_t j = 0; j < coefs[i].size(); ++j)
			secondary_dist.push_back(unify_legandre_polinom(coefs[i][j]));
		legandre_coefs.push_back(std::move(secondary_dist));
	}
}

std::vector<double> sum_vectors_with_coefs(
	const std::vector<double>& vec1, double coef1,
	const std::vector<double>& vec2, double coef2) {
	std::vector<double> result(vec1.size());
	for (size_t k = 0; k < result.size(); ++k)
		result[k] = vec1[k] * coef1 + vec2[k] * coef2;
	return result;
}

std::vector<double> sum_vectors_with_coefs(
	std::vector<double>&& vec1, double coef1,
	std::vector<double>&& vec2, double coef2) {
	std::vector<double> result(std::move(vec1));
	for (size_t k = 0; k < result.size(); ++k)
		result[k] = result[k] * coef1 + vec2[k] * coef2;
	return result;
}

std::vector<double> Distribution2D::interpolate(int i, double point) const {
	int j = static_cast<int>(point);
	if (i < 0 || i >= legandre_coefs.size())
		return std::vector<double>(dim);
	if (j < 0 || j >= legandre_coefs[i].size() - 1)
		return std::vector<double>(dim);

	double ratio = point - j;
	return sum_vectors_with_coefs(legandre_coefs[i][j], 1 - ratio, legandre_coefs[i][j + 1], ratio);
}

std::vector<double> Distribution2D::interpolate(double primary_point, double point1, double point2) const {
	int i = static_cast<int>(primary_point);
	double ratio = primary_point - i;

	return sum_vectors_with_coefs(
		interpolate(i, point1), 1 - ratio, 
		interpolate(i + 1, point2), ratio);
}

#endif
