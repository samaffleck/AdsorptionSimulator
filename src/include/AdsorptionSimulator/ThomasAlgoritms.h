#pragma once

#include "Eigen/Dense"

#include <vector>

namespace LinearSolver
{
	void TDMA(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V);
}


class ThomasAlgorithm {

public:

	ThomasAlgorithm() = default;
	~ThomasAlgorithm() = default;

	int n{};
	std::vector<double> cPrime{};
	std::vector<double> dPrime{};

	void initialise(int size);
	void solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& X, Eigen::VectorXd& V);
	void solve(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V);
	void printVector(const Eigen::VectorXd& v) const;

	static double calculateError(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& Sp, const Eigen::VectorXd& V);

};
