#include "Experiment.hpp"

#ifdef WIN32
std::optional<Experiment> Experiment::load_from_dll(std::filesystem::path path) noexcept {
	return std::nullopt;
}
#else
std::optional<Experiment> Experiment::load_from_dll(std::filesystem::path path) noexcept {
	return std::nullopt;
}
#endif

