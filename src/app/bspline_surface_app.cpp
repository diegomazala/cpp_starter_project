#include <iostream>
#include <array>
#include "bspline_surface.h"


int main(int argc, char* argv[])
{
	using decimal_t = double;
	using vec3_t = std::array<decimal_t, 3>;

	const size_t width = 16;
	const size_t height = width;

	const uint32_t m = (argc > 1 ? atoi(argv[1]) : 64);
	const uint32_t n = m;

	std::vector<vec3_t> points(width * height);
	for (size_t i = 0; i < height; ++i)
	{
		for (size_t j = 0; j < width; ++j)
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

	std::cout << std::fixed  << "Error: " << sum_error << std::endl;


    return EXIT_SUCCESS;
}

