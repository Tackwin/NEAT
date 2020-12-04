#pragma once

#include <optional>
#include <vector>
#include <array>
#include "AI/Neat.hpp"

struct Network_Window {
	std::vector<float> inputs;
	std::vector<float> outputs;

	std::vector<size_t> links_created;
	std::vector<size_t> links_deleted;
	std::vector<size_t> links_hovered;

	std::optional<float> new_weight;

	void embed_render(Network& network) noexcept;
};

struct Genome_Window {
	void embed_render(Genome& genome) noexcept;
};

struct Neat_Window {
	bool open_edit_initial_network = false;
	bool open_best_network = false;

	bool run_generation = false;
	bool auto_run = false;

	Network_Window initial_network_window;
	Genome_Window initial_genome_window;
	Network initial_network;
	Genome initial_genome;

	Network_Window best_network_window;

	enum Evaluation {
		XOR
	} evaluation = Evaluation::XOR;

	std::vector<float> max_fitness;
	std::vector<float> max_adjusted_fitness;
	std::vector<std::array<float, 100>> fitness_histograms;
	std::vector<std::array<float, 100>> adjusted_fitness_histograms;
	size_t max_fitness_n_samples = 100;

	void render(Neat& neat) noexcept;

	void get_stats(const std::vector<Neat::Result>& results) noexcept;
};

void update_genome(Genome& gen, const Network_Window& network_window_input) noexcept;