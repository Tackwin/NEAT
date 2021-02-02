#pragma once

#include <optional>
#include <vector>
#include <array>
#include "AI/Neat.hpp"
#include "imgui/imnodes.h"

// Like wtf imnodes ??? why don't you define your struct in your header??
// I'm now forced to explicitely implement the full combination of move, copy constructor
// destructor ?????
struct EditorContext {
	imnodes::EditorContext* ctx;

	EditorContext() { ctx = imnodes::EditorContextCreate(); }
	EditorContext(EditorContext&&) { ctx = imnodes::EditorContextCreate(); }
	~EditorContext() { imnodes::EditorContextFree(ctx); }
};

struct Network_Window {
	EditorContext ctx;

	std::vector<float> inputs;
	std::vector<float> outputs;

	std::vector<size_t> links_created;
	std::vector<size_t> links_deleted;
	std::vector<size_t> links_hovered;

	std::optional<float> new_weight;

	bool hovered = false;
	bool layout = false;

	void embed_render(Network& network) noexcept;
	void auto_layout(Network& network) noexcept;
};

struct Genome_Window {
	void embed_render(Genome& genome) noexcept;
};

struct Population_Window {
	size_t generation = 0;
	float min_fitness = 0;
	float max_fitness = 0;
	size_t min_neurons = 0;
	size_t max_neurons = 0;

	std::vector<Network_Window> network_windows;

	void render(
		const std::vector<std::vector<Genome>>& pop_snapshots,
		const std::vector<std::vector<Neat::Result>>& res_snapshots
	) noexcept;
	void embed_render(
		const std::vector<std::vector<Genome>>& pop_snapshots,
		const std::vector<std::vector<Neat::Result>>& res_snapshots
	) noexcept;

};

struct Neat_Window {
	bool open_edit_initial_network = false;
	bool open_best_phenotype = false;
	bool open_best_network = false;

	bool run_generation = false;
	bool auto_run = false;

	Population_Window population_window;

	Network_Window initial_network_window;
	Genome_Window initial_genome_window;
	Network initial_network;
	Genome initial_genome;

	Network_Window best_network_window;
	Genome_Window best_genome_window;
	Genome best_genome;

	void (*render_phenotype)(Network&, void*);

	enum Evaluation {
		XOR,
		SPV,  // Single pole velocity
		DPV,  // Double pole velocity
		DP,   // Double pole no velocity
		HDPV  // Double pole velocity
	} evaluation = Evaluation::DPV;

	std::vector<float> max_fitness;
	std::vector<float> max_adjusted_fitness;
	std::vector<std::array<float, 100>> fitness_histograms;
	std::vector<std::array<float, 100>> adjusted_fitness_histograms;
	size_t max_fitness_n_samples = 100;

	struct Specie_Info {
		size_t n  = 0;
		size_t id = 0;
	};
	std::vector<std::vector<Specie_Info>> species_infos;

	bool capture_population = true;
	std::vector<std::vector<Genome>> population_snapshots;
	std::vector<std::vector<Neat::Result>> results_snapshots;

	bool open_explore_population = false;

	void render(Neat& neat) noexcept;
	void render_best(Neat& neat) noexcept;

	void get_stats(const Neat& neat, const std::vector<Neat::Result>& results) noexcept;
};


void update_genome(Genome& gen, const Network_Window& network_window_input) noexcept;