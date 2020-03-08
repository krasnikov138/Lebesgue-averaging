#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <list>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

struct Segment {
	double l, r;
	Segment(double l = 0, double r = 0): l(l), r(r) {};
};

typedef std::vector<Segment> Carrier;

struct Spectrum {
public:
	std::string name;
	size_t size;
	std::vector<double> e, cs;

	Spectrum(): name("Unnamed"), size(0) {};
	Spectrum(std::string name): name(name), size(0) {};

	void load(std::string file_name);

	double get_energy(double point);
	double get_point(double energy, int index);
	double get_cs(double point);
	double get_cs(int i, double ratio);
	double linear(double energy, int i);
	double resonance_begin(double threshold, int use_points);
};

void Spectrum::load(std::string file_name) {
	std::ifstream input(file_name);
	double e_val, cs_val;
	while (input) {
		input >> e_val >> cs_val;
		e.push_back(e_val);
		cs.push_back(cs_val);
	}
	size = e.size();
}

double Spectrum::get_energy(double point) {
	int i = static_cast<int>(point);
	double ratio = point - i;
	return ratio * e[i + 1] + (1 - ratio) * e[i];
}

double Spectrum::get_point(double energy, int index) {
	double point = index;
	if (index + 1 < size)
		point += (energy - e[index]) / (e[index + 1] - e[index]);
	return point;
}

double Spectrum::get_cs(int i, double ratio) {
	return ratio * cs[i + 1] + (1 - ratio) * cs[i];
}

double Spectrum::get_cs(double point) {
	int i = static_cast<int>(point);
	double ratio = point - i;
	return get_cs(i, ratio);
}

double Spectrum::linear(double energy, int i) {
	double ratio = (energy - e[i]) / (e[i + 1] - e[i]);
	return get_cs(i, ratio);
}

double Spectrum::resonance_begin(double threshold = 5.0, int use_points = 10) {
	double der_mean = 0;
	size_t i = 1;
	while (i < size) {
		double der = fabs((log10(cs[i]) - log10(cs[i - 1])) / (log10(e[i]) - log10(e[i - 1])));
		if (i < use_points) 
			der_mean += der / use_points;
		else if (der > threshold * der_mean)
			break;
		i++;
	}
	return e[i];
}

bool is_close(double a, double b, double rel_tol=1e-8) {
	if (fabs(a - b) < rel_tol * b) 
		return true;
	return false;
}

bool is_less(double a, double b, double rel_tol=1e-8) {
	if (a - b < rel_tol * b)
		return true;
	return false;
}

int find_dominant_mat(std::vector<Spectrum> &spectra, const std::vector<int> &offset, double e) {
	int index_max = -1;
	double cs_max = 0.0;
	for (int i = 0; i < spectra.size(); i++) {
		double cs = spectra[i].linear(e, offset[i]);
		if (cs_max < cs) {
			index_max = i;
			cs_max = cs;
		}
	}
	return index_max;
}

std::vector<std::vector<Carrier>> build_carriers(std::vector<Spectrum> &spectra, 
	double start_energy = 1e-5, double K = 0.1) {

	std::vector<int> offset(spectra.size(), 0);
	std::vector<std::list<Segment>> segments(spectra.size());
	int mat, prev_mat;

	double e = start_energy;
	/* skip unnecessary energy points */
	for (int i = 0; i < spectra.size(); i++) {
		while (is_less(spectra[i].e[offset[i]], start_energy))
			offset[i]++;
		offset[i]--;
	}

	double pos, prev_pos;
	/* find first dominant material and relative spectrum position */
	prev_mat = find_dominant_mat(spectra, offset, e);
	prev_pos = spectra[mat].get_point(e, offset[prev_mat]);

	for (auto &p: offset) p++;
	while (e > 0) {
		/* next energy point */
		e = -1.0;
		for (int i = 0; i < spectra.size(); i++) {
			if (offset[i] < spectra[i].size) {
				double energy = spectra[i].e[offset[i]];
				if (e < 0 || e > energy)
					e = energy;
			}
		}
		if (e > 0) {
			mat = find_dominant_mat(spectra, offset, e);
			pos = spectra[mat].get_point(e, offset[mat] - 1);
			if (mat == prev_mat) {
				if (segments[mat].size() > 0 && is_close(floor(pos), floor(segments[mat].back().l)))
					segments[mat].back().r = pos;
				else
					segments[mat].push_back(Segment(prev_pos, pos));
			}
			prev_mat = mat;
			prev_pos = pos;
		}
		/* skip close points */
		for (int i = 0; i < spectra.size(); i++) {
			if (offset[i] < spectra[i].size && is_close(e, spectra[i].e[offset[i]]))
				offset[i]++;
		}
	}

	std::vector<std::vector<Carrier>> res;
	for (int i = 0; i < spectra.size(); i++) {
		Carrier carrier;
		std::vector<Carrier> parts;
		double begin;
		while (!segments[i].empty()) {
			Segment current = segments[i].front();
			if (carrier.empty()) 
				begin = current.l;
			if ((current.r - begin) / (current.r + begin) < K / 2) {
				carrier.push_back(current);
				segments[i].pop_front();
			} else if ((current.l - begin) / (current.l + begin) > K / 2) {
				parts.push_back(carrier);
				carrier.clear();
			} else {
				double end = (1 + K / 2) / (1 - K / 2) * begin;
				carrier.push_back(Segment(current.l, end));
				segments[i].front().l = end;
				parts.push_back(carrier);
				carrier.clear();
			}
		}
		res.push_back(parts);
	}

	return res;
}

double measure(double s, Spectrum& spectrum, Carrier& carrier) {
	double mes = 0;
	double e_l, e_r, cs_l, cs_r;

	for (int i = 0; i < carrier.size(); i++) {
		e_l = spectrum.get_energy(carrier[i].l);
		e_r = spectrum.get_energy(carrier[i].r);
		cs_l = spectrum.get_cs(carrier[i].l);
		cs_r = spectrum.get_cs(carrier[i].r);

		if (cs_l < s && cs_r < s)
			mes += (e_r - e_l);
		else if (cs_l < s && cs_r >= s)
			mes += (e_r - e_l) * (s - cs_l) / (cs_r - cs_l);
		else if (cs_l >= s && cs_r < s)
			mes += (e_r - e_l) * (s - cs_r) / (cs_l - cs_r);
	}
	return mes;
}

std::tuple<double, double> cross_section_interval(Spectrum& spectrum, Carrier& carrier) {
	double cs_min = 1.0e9;
	double cs_max = 0.0;
	double cs_l, cs_r;
	for (int i = 0; i < carrier.size(); i++) {
		cs_l = spectrum.get_cs(carrier[i].l);
		cs_r = spectrum.get_cs(carrier[i].r);
		cs_min = std::min(cs_min, std::min(cs_l, cs_r));
		cs_max = std::max(cs_max, std::max(cs_l, cs_r));
	}
	return std::make_tuple(cs_min, cs_max);
}

int main(int argc, char** argv) {
	/* read command lines arguments */
	double K = 0.1;
	if (argc > 2) {
		std::cout << "Usage: ./build_resonances [K]";
		exit(0);
	} else if (argc == 2)
		K = std::stod(argv[1]);

	/* load all files from material directory */
	std::vector<Spectrum> spectra;
	const std::string mat_directory = "materials";
	for (auto& entry : fs::directory_iterator(mat_directory)) {
		std::string name(entry.path().filename());

		Spectrum sp(name);
		sp.load(entry.path());
		spectra.push_back(sp);
	}

	double start_energy = 1e8;
	for (auto &s : spectra) {
		double start = s.resonance_begin();
		if (start < start_energy)
			start_energy = start;
	}
	std::vector<std::vector<Carrier>> res = build_carriers(spectra, start_energy, K);

	fs::path dir_name("carriers");
	/* write carriers in files */
	fs::remove_all(dir_name);
	for (int i = 0; i < res.size(); i++) {
		fs::create_directories(dir_name / spectra[i].name);
		for (int j = 0; j < res[i].size(); j++) {
			std::string fname = std::to_string(j) + ".csv";
			std::ofstream f(dir_name / spectra[i].name / fname);
			f << std::setprecision(10);
			for (int k = 0; k < res[i][j].size(); k++) {
				f << spectra[i].get_energy(res[i][j][k].l) << ',';
				f << spectra[i].get_cs(res[i][j][k].l) << ',';
				f << spectra[i].get_energy(res[i][j][k].r) << ',';
				f << spectra[i].get_cs(res[i][j][k].r) << '\n';
			}
		}
	}

	/* calculate measure variable for some cases */
	const int N = 100;
	dir_name = "measures";
	fs::remove_all(dir_name);
	for (int i = 0; i < res.size(); i++) {
		fs::create_directories(dir_name / spectra[i].name);
		for (int j = 0; j < res[i].size(); j++) {
			std::string fname = std::to_string(j) + ".csv";
			std::ofstream f(dir_name / spectra[i].name / fname);
			f << std::setprecision(10);
			auto [cs_min, cs_max] = cross_section_interval(spectra[i], res[i][j]);
			for (int k = 0; k < N + 1; k++) {
				double s = cs_min + (cs_max - cs_min) * k / N;
				f << s << ',' << measure(s, spectra[i], res[i][j]) << '\n';
			}
		}
	}

	return 0;
}
