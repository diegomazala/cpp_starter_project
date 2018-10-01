#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

// bspline namespace starts here
namespace bspline
{

	//
	// B-splines functions
	//
	template <typename decimal_t>
	static std::function<decimal_t(decimal_t)> B[4] = {
		[](decimal_t t) { return std::pow((1 - t), 3) / 6; },
		[](decimal_t t) {
			return (3 * std::pow(t, 3) - 6 * std::pow(t, 2) + 4) / 6;
		},
		[](decimal_t t) {
			return (-3 * std::pow(t, 3) + 3 * std::pow(t, 2) + 3 * t + 1) / 6;
		},
		[](decimal_t t) { return std::pow(t, 3) / 6; } };

	//
	// Function to test doubles and floats considering decimals
	//
	template <typename T>
	constexpr bool logically_equal(const T &a, const T &b,
		const T error_factor = 1)
	{
		return a == b ||
			std::abs(a - b) < std::abs(std::min(a, b)) *
			std::numeric_limits<T>::epsilon() * error_factor;
	}

	template <typename decimal_t, typename vec3_t = std::array<decimal_t, 3>>
	class surface
	{
	public:

		using matrix_t = std::vector<std::vector<decimal_t>>;
		using matrix4_t = std::array<std::array<decimal_t, 4>, 4>;
		// using vec3_t = std::array<decimal_t, 3>;
		using vec2_t = std::array<decimal_t, 2>;

		surface(const vec3_t* _points, size_t _point_count, uint32_t _m, uint32_t _n)
			: points(_points), point_count(_point_count), m(_m), n(_n), z_surface(point_count)
		{
			this->init();
			this->apply_offset();
			this->compute();
			this->remove_offset();
		}

		surface() = delete;

		//
		// Run the algorith to build bspline surface
		//
		void compute()
		{
			matrix_t delta = { (m + 3), std::vector<decimal_t>(n + 3, 0) };
			matrix_t omega = { (m + 3), std::vector<decimal_t>(n + 3, 0) };
			this->phi = { (m + 3), std::vector<decimal_t>(n + 3, 0) };

			const decimal_t interval_normalization_factor_u = m * urange_inv;
			const decimal_t interval_normalization_factor_v = n * vrange_inv;

			auto points_ptr = points;
			for (size_t index = 0; index < point_count; ++index, points_ptr++)
			{
				const decimal_t x = (*points_ptr)[0];
				const decimal_t y = (*points_ptr)[1];
				const decimal_t z = z_surface[index]; //p[2];   --> z_surface is offseted with average_z

				// Map to the half open domain Omega = [0,m) x [0,n)
				// The mapped u and v must be (strictly) less than m and n respectively
				decimal_t u = (x - umin) * interval_normalization_factor_u;
				decimal_t v = (y - vmin) * interval_normalization_factor_v;

				//
				// Compute [i,j] indices and [s,t] params
				//
				auto[i, j, s, t] = compute_ijst(u, v);
				i += 1;	// +1 to go for positive indices
				j += 1;

				//
				// Compute w_kl's and sum_sum w²_ab
				//
				const auto[w, sum_w_ab2] = compute_wkl_and_sum_w_ab2(s, t);
				const auto sum_w_ab2_inv = 1 / sum_w_ab2;

				//
				// Update delta and gamma matrices
				//
				for (auto k = 0; k < 4; ++k)
				{
					for (auto l = 0; l < 4; ++l)
					{
						// eq(3)
						auto phi_kl = w[k][l] * z * sum_w_ab2_inv;

						auto w2 = w[k][l] * w[k][l];

						delta[i + k][j + l] += w2 * phi_kl;
						omega[i + k][j + l] += w2;
					}
				}
			}

			for (uint32_t i = 0; i < m + 3; ++i)
			{
				for (uint32_t j = 0; j < n + 3; ++j)
				{
					if (!logically_equal(omega[i][j], static_cast<decimal_t>(0)))
					{
						this->phi[i][j] = delta[i][j] / omega[i][j];
					}
					else
					{
						this->phi[i][j] = 0;
					}
				}
			}
		}

		//
		// Evaluate the function at (u, v)
		//
		decimal_t eval(decimal_t u, decimal_t v) const 
		{
			// Map to the half open domain Omega = [0,m) x [0,n)
			u = (u - umin) / (umax - umin) * m;
			v = (v - vmin) / (vmax - vmin) * n;

			//
			// Compute [i,j] indices and [s,t] params
			//
			auto[i, j, s, t] = compute_ijst(u, v);

			i += 1;    // +1 to go for positive indices
			j += 1;

			//
			// Evaluate the function
			//
			decimal_t val = 0;
			for (auto k = 0; k < 4; ++k)
			{
				for (auto l = 0; l < 4; ++l)
				{
					val +=
						this->phi[i + k][j + l] * B<decimal_t>[k](s) * B<decimal_t>[l](t);
				}
			}
			return val;
		}

		//
		// Overload operator () to evaluate the function at (u, v)
		//
		decimal_t operator () (decimal_t u, decimal_t v) const
		{
			return eval(u, v);
		}

		//
		// Compute average for z and use it for initialize the surface
		//
		void apply_offset()
		{
			auto points_ptr = points;
			for (size_t index = 0; index < point_count; ++index, points_ptr++)
			{
				z_surface[index] = (*points_ptr)[2] - average_z;
			}
		}

		// 
		// Remove the initial offset applied to initialize the surface
		//
		void remove_offset()
		{
			for (auto& v1 : this->phi)
			{
				for (auto& v2 : v1)
				{
					v2 += average_z;
				}
			}
			//average_z = 0;
		}


	protected:

		//
		// Compute [i,j] indices and [s,t] params
		//
		std::tuple<int64_t, int64_t, decimal_t, decimal_t> compute_ijst(decimal_t x,
			decimal_t y) const
		{
			vec2_t floor = { std::floor(x), std::floor(y) };
			auto i = static_cast<int64_t>(floor[0] - 1);
			auto j = static_cast<int64_t>(floor[1] - 1);
			decimal_t s = x - floor[0];
			decimal_t t = y - floor[1];

			// test if with uv bigger (0.1), this is still needed
			if (i == m - 1)
			{
				i--;
				s = 1;
			}
			if (j == n - 1)
			{
				j--;
				t = 1;
			}

			return { i, j, s, t };
		}

		//
		// Compute w_kl's and sum_sum w²_ab
		//
		static std::tuple<matrix4_t, decimal_t>
			compute_wkl_and_sum_w_ab2(decimal_t s, decimal_t t)
		{
			matrix4_t w;
			decimal_t sum_w_ab2 = 0;
			for (auto k = 0; k < 4; ++k)
			{
				for (auto l = 0; l < 4; ++l)
				{
					w[k][l] = B<decimal_t>[k](s) * B<decimal_t>[l](t);
					sum_w_ab2 += w[k][l] * w[k][l];
				}
			}
			return { w, sum_w_ab2 };
			// return std::make_tuple(w, sum_w_ab2);
		}




		void init()
		{
			umin = std::numeric_limits<decimal_t>::max();
			vmin = std::numeric_limits<decimal_t>::max();
			umax = std::numeric_limits<decimal_t>::min();
			vmax = std::numeric_limits<decimal_t>::min();

			decimal_t sum_z = 0;

			auto points_ptr = points;
			for (size_t index = 0; index < point_count; ++index, points_ptr++)
			{
				const decimal_t x = (*points_ptr)[0];
				const decimal_t y = (*points_ptr)[1];
				const decimal_t z = (*points_ptr)[2];

				if (x < umin) { umin = x; }
				if (y < vmin) { vmin = y; }
				if (x > umax) { umax = x; }
				if (y > vmax) { vmax = y; }

				sum_z += z;
			}

			urange_inv = 1 / (umax - umin);
			vrange_inv = 1 / (vmax - vmin);

			average_z = sum_z / point_count;
		}

		//
		// Attributes
		const vec3_t* points;
		const size_t point_count;
		uint32_t m, n;
		matrix_t phi;
		decimal_t umin, vmin, umax, vmax, urange_inv, vrange_inv, average_z;
		std::vector<decimal_t> z_surface;
	};

} // namespace bspline
