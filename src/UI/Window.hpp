#pragma once

#include "AI/Neat.hpp"

struct Network_Window {
	std::vector<float> inputs;
	std::vector<float> outputs;

	std::vector<size_t> links_created;
	std::vector<size_t> links_deleted;

	void embed_render(Network& network) noexcept;
};

struct Genome_Window {
	void embed_render(Genome& genome) noexcept;
};

struct Neat_Window {

	size_t population = 0;

	void render(Neat& neat) noexcept;
};

void update_genome(Genome& gen, const Network_Window& network_window_input) noexcept;