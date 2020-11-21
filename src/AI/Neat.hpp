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
	struct Node {
		float y;
		enum Func {
			Relu
		} f;
	};

	std::vector<Connections> connections;
	std::vector<Node> nodes;
	std::vector<size_t> input_nodes;
	std::vector<size_t> output_nodes;

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

	std::vector<Node_Gene> node_genes;
	std::vector<Connect_Gene> connect_genes;

	static Genome mate(const Genome& a, const Genome& b) noexcept;

	void add_connection(std::uint32_t in, std::uint32_t out, float w) noexcept;
	void add_new_input() noexcept;
	void add_new_output() noexcept;
	void add_node() noexcept;

	void del_connection(size_t i) noexcept;

	Network phenotype() noexcept;
};

struct Neat {
	size_t n_inputs;
	size_t n_outputs;
	std::vector<Genome> population;

	void fill_with(Genome g) noexcept;

	void add_agents(size_t n) noexcept;

	void evaluate(std::function<float(Network&)> fitness) noexcept;
	void turnament_evaluate(
		std::function<float(Network&)> fitness, size_t game_size
	) noexcept;
};