#pragma once
#include <vector>
#include <atomic>
#include <functional>
#include <cstdint>

static std::atomic<std::uint32_t> Innov_N = 0;

enum Node_Activation_Func {
	Relu = 0,
	Lin,
	Sig,
	Count
};

static const char* node_activation_func_str(Node_Activation_Func x) noexcept {
	#define X(x) case x: return #x;
	switch (x) {
	X(Relu);
	X(Lin);
	X(Sig);
	default:
		return "?";
	}
	#undef X
}

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
		float y = 0;
		float s = 0;
		Node_Activation_Func f = Node_Activation_Func::Relu;
	};

	std::vector<Connections> connections;
	std::vector<Node> nodes;
	std::vector<Node_Kind> node_kinds;

	void compute(float* in, size_t in_size, float* out, size_t out_size) noexcept;

};

struct Genome {
	struct Node_Gene {
		std::uint32_t node_id = 0;
		Node_Activation_Func f = Node_Activation_Func::Relu;
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
	static float dist(
		const Genome& a, const Genome& b, float c1 = 1, float c2 = 1, float c3 = 1
	) noexcept;

	void add_connection(std::uint32_t in, std::uint32_t out, float w) noexcept;
	void add_new_input() noexcept;
	void add_new_hidden() noexcept;
	void add_new_output() noexcept;
	void add_node(size_t c) noexcept;
	void del_node(size_t c) noexcept;

	void del_connection(size_t i) noexcept;
	void dis_connection(size_t i) noexcept;
	void update_weight(size_t i, float w) noexcept;

	Network phenotype() noexcept;

	std::vector<std::uint8_t> serialize() const noexcept;
	static Genome deserialize(const std::vector<std::uint8_t>& data) noexcept;
};

struct Neat {
	float young_advantage = 1.f;
	float survival_rate = 0.5f;
	float mutation_rate = 0.1f;
	float specie_crowd_rate = 1.f;
	float complexity_cost = 0.1f;

	float c_1 = 1.0f;
	float c_2 = 1.0f;
	float c_3 = 0.4f;
	float specie_dt = 3.f;
	size_t preferred_number_of_species = 0;

	bool specifie_number_of_species = false;

	size_t generation_number = 0;
	
	size_t population_size = 10'000;
	std::vector<Genome> population;

	struct Genome_Info {
		size_t age;
		size_t specie;
	};
	std::vector<Genome_Info> genome_info;

	struct Specie {
		inline static size_t Specie_N = 0;

		Genome repr;
		size_t size = 0;
		size_t num = 0;
		size_t gen_since_upgrade = 0;

		float best_fitness = FLT_MIN;

		std::vector<size_t> idx_in_population;
	};
	std::vector<Specie> species;


	struct Result {
		Genome* g;
		float fitness = 0;
		float adjusted_fitness = 0;
	};
	std::vector<Result> results;

	void complete_with(Genome g, size_t to) noexcept;

	void add_genome(Genome g) noexcept;

	void evaluate(std::function<float(Network&)> fitness) noexcept;
	void speciate() noexcept;
	void select() noexcept;
	void populate() noexcept;

	Genome mutation(const Genome& in) noexcept;

	void turnament_evaluate(

		std::function<float(Network&)> fitness, size_t game_size
	) noexcept;
};

extern float xor_fitness(Network& net) noexcept;