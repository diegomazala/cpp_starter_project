#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include <array>
#include "bspline_surface.h"

using decimal_t = double;

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


decimal_t base_test()
{
	using vec3_t = std::array<decimal_t, 3>;

	const size_t width = 16;
	const size_t height = 16;

	const uint32_t m = 64;
	const uint32_t n = 64;

	std::vector<vec3_t> points(width * height);
	for (auto i = 0; i < height; ++i)
	{
		for (auto j = 0; j < width; ++j)
		{
			auto index = i * width + j;
			points[index][0] = static_cast<decimal_t>(i) / static_cast<decimal_t>(width);
			points[index][1] = static_cast<decimal_t>(j) / static_cast<decimal_t>(height);
			points[index][2] = std::sin(points[index][0] * 2);
		}
	}

	bspline::surface<decimal_t, vec3_t> surf = { points.data(), points.size(), m, n };

	decimal_t sum_error = 0;
	for (auto p : points)
	{
		auto f = surf(p[0], p[1]);
		sum_error += std::abs(f - p[2]);
	}

	return sum_error;
}




TEST_CASE("BSplineSurface", "[bspline_surface]") 
{
	REQUIRE(logically_equal<decimal_t>(base_test(), 0));
}




