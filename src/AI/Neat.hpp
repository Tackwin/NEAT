#pragma once
#include <vector>
#include <atomic>
#include <functional>
#include <cstdint>

static std::atomic<std::uint32_t> Innov_N = 0;


struct Network {
	struct Connections {
		std::uint32_t in;
		std::uint32_t out;
		float w;
	};
	enum Node_Kind {
		Sensor,
		Output,
		Hidden
	};
	struct Node {
		float y;
		enum Func {
			Relu
		} f;
	};

	std::vector<Connections> connections;
	std::vector<Node> nodes;
	std::vector<Node_Kind> node_kinds;

	void compute(float* in, size_t in_size, float* out, size_t out_size) noexcept;
};

struct Genome {
	struct Node_Gene {
		std::uint32_t node_id = 0;
		enum Kind {
			Sensor,
			Output,
			Hidden
		} kind = Kind::Sensor;
	};

	struct Connect_Gene {
		std::uint32_t in_id  = 0;
		std::uint32_t out_id = 0;
		std::uint32_t innov  = 0;
		float w = 0;
		bool enabled = true;
	};

	size_t input_nodes = 0;
	size_t output_nodes = 0;

	std::vector<Node_Gene> node_genes;
	std::vector<Connect_Gene> connect_genes;

	static Genome mate(const Genome& a, const Genome& b) noexcept;

	void add_connection(std::uint32_t in, std::uint32_t out, float w) noexcept;
	void add_new_input() noexcept;
	void add_new_hidden() noexcept;
	void add_new_output() noexcept;
	void add_node(size_t c) noexcept;

	void del_connection(size_t i) noexcept;
	void dis_connection(size_t i) noexcept;
	void update_weight(size_t i, float w) noexcept;

	Network phenotype() noexcept;
};

struct Neat {
	size_t n_inputs = 0;
	size_t n_outputs = 0;

	float survival_rate = 0.5f;

	float mutation_rate = 0.1f;

	size_t population_size = 0;
	std::vector<Genome> population;

	struct Result {
		Genome* g;
		float fitness;
	};
	std::vector<Result> results;

	void fill_with(Genome g) noexcept;

	void add_agents(size_t n) noexcept;

	void evaluate(std::function<float(Network&)> fitness) noexcept;
	void select() noexcept;
	void populate() noexcept;

	void turnament_evaluate(
		std::function<float(Network&)> fitness, size_t game_size
	) noexcept;
};