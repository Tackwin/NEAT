#pragma once

#include <filesystem>
#include <optional>
#include "AI/Neat.hpp"

struct Experiment {
	using Ftiness_f = float(*)(Network&, void*);
	using Render_f = void(*)(Network&, void*);
	using Alloc_f = void*(*)(Network&);
	using Free_f = void(*)(Network&, void*);

	using Genome_f = Genome(*)();

	std::string name = "";

	Genome_f genome = nullptr;

	Ftiness_f fitness = nullptr;
	Render_f render = nullptr;
	Alloc_f alloc = nullptr;
	Free_f free = nullptr;



	static std::optional<Experiment> load_from_dll(std::filesystem::path path) noexcept;
};