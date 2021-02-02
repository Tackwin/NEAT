#pragma once

#include <chrono>

namespace xstd {

	static double seconds() noexcept {
		auto now     = std::chrono::system_clock::now();
		auto epoch   = now.time_since_epoch();
		auto seconds = std::chrono::duration_cast<std::chrono::nanoseconds>(epoch);

		// return the number of seconds
		return seconds.count() / 1'000'000'000.0;
	}

};
static double randomd() noexcept { return (double)rand() / (1.0 + RAND_MAX); }
static double randomd2() noexcept { auto d = randomd(); return d * d; }

static inline float sig(float x)  noexcept { return -1 + 2 / (1 + std::expf(-1.f * x)); }
static inline float per(float x)  noexcept { return x > 0.5 ? 1 : 0; }
static inline float relu(float x) noexcept { return x > 0 ? x : 0; }
static inline float lin(float x)  noexcept { return x; }