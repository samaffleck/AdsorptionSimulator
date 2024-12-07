#include "AdsorptionSimulator/ThomasAlgoritms.h"


namespace LinearSolver
{
	void TDMA(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V)
	{
		// Forward pass
		auto n = int(V.size());
		std::vector<double> cPrime(n);
		std::vector<double> dPrime(n);

		cPrime[0] = Ae(0) / Ap(0);
		dPrime[0] = S[0] / Ap(0);
		for (int i = 1; i < n - 1; i++)
		{
			cPrime[i] = Ae(i) / (Ap(i) - Aw(i) * cPrime[i - 1]);
			dPrime[i] = (S[i] - Aw(i) * dPrime[i - 1]) / (Ap(i) - Aw(i) * cPrime[i - 1]);
		}
		dPrime[n - 1] = (S[n - 1] - Aw(n - 1) * dPrime[n - 2]) / (Ap(n - 1) - Aw(n - 1) * cPrime[n - 2]);

		// Back substitution
		V[n - 1] = dPrime[n - 1];
		for (int i = n - 2; i >= 0; i--)
		{
			V[i] = dPrime[i] - cPrime[i] * V[i + 1];
		}
	}

	double getError(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V)
	{
		const auto N = (int)V.size();
		Eigen::VectorXd residual;
		residual.setConstant(N, 0);

		residual(0) = Ap(0) * V(0) + Ae(0) * V(1) - S(0);

		for (size_t n = 1; n < N - 1; ++n)
		{
			residual(n) = Aw(n) * V(n - 1) + Ap(n) * V(n) + Ae(n) * V(n + 1) - S(n);
		}

		residual(N - 1) = Aw(N - 1) * V(N - 2) + Ap(N - 1) * V(N - 1) - S(N - 1);

		return residual.norm();
	}
}
