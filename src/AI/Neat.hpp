#pragma once
#include <unordered_map>
#include <vector>
#include <atomic>
#include <functional>
#include <optional>
#include <cstdint>
#include <float.h>

static std::atomic<std::uint32_t> Innov_N = 0;
static std::atomic<size_t> Genome_Serial = 0;

enum Node_Activation_Func {
	Sig = 0,
	Per,
	Lin,
	Relu,
	Count
};

static const char* node_activation_func_str(Node_Activation_Func x) noexcept {
	#define X(x) case x: return #x;
	switch (x) {
	X(Relu);
	X(Per);
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
		Node_Activation_Func f = Node_Activation_Func::Sig;
	};

	std::vector<Connections> connections;
	std::vector<Node> nodes;
	std::vector<Node_Kind> node_kinds;

	size_t input_nodes = 0;
	size_t output_nodes = 0;

	size_t genome_serial_number = 0;

	void clear() noexcept;
	void compute(float* in, size_t in_size, float* out, size_t out_size) noexcept;
	void compute_clear(float* in, size_t in_size, float* out, size_t out_size) noexcept;
};

struct Genome {
	struct Node_Gene {
		std::uint32_t node_id = 0;
		Node_Activation_Func f = Node_Activation_Func::Sig;
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

	size_t serial_number = 0;

	static Genome mate(const Genome& a, const Genome& b) noexcept;
	static bool dist(
		const Genome& a,
		const Genome& b,
		float c1 = 1,
		float c2 = 1,
		float c3 = 1,
		float tresh = 3
	) noexcept;

	void add_connection(std::uint32_t in, std::uint32_t out, float w) noexcept;
	void add_connection(const Connect_Gene& other) noexcept;
	void add_new_input() noexcept;
	void add_new_hidden() noexcept;
	void add_new_output() noexcept;
	void add_node(size_t c) noexcept;
	void del_node(size_t c) noexcept;

	void del_connection(size_t i) noexcept;
	void dis_connection(size_t i) noexcept;
	void update_weight(size_t i, float w) noexcept;

	Network phenotype() const noexcept;

	std::vector<std::uint8_t> serialize() const noexcept;
	static Genome deserialize(const std::vector<std::uint8_t>& data) noexcept;
};

struct Neat {
	struct Mutation_Params {
		float new_connection_weight_rate = 0.1f;
		float update_connection_rate     = 0.8f;
		float toggle_connection_rate    = 0.1f;
		float add_connection_rate        = 0.3f;
		float update_node_act_rate       = 0.1f;
		float add_node_rate              = 0.01f;
		float del_node_rate              = 0.01f;
	} params;

	float young_advantage = 0.f;
	float survival_rate = 0.5f;
	float mutation_rate = 1.0f;
	float specie_crowd_rate = 1.f;
	float complexity_cost = 0.5f;
	float population_competition_rate = 0.0f;

	size_t age_cutoff_sterile_specie = 20;
	size_t min_specie_size_advantage = 100;

	float c_1 = 1.0f;
	float c_2 = 1.0f;
	float c_3 = 0.4f;
	float specie_dt = 3.f;
	size_t preferred_number_of_species = 20;

	bool specifie_number_of_species = false;

	size_t generation_number = 0;
	
	size_t population_size = 10'000;
	std::vector<Genome> population;

	// >SEE(Tackwin): this doesn't exactly concern Neat where should i put it ??
	size_t threads_to_use = 4;

	struct Genome_Info {
		size_t age;
		std::optional<size_t> specie;
	};
	std::vector<Genome_Info> genome_info;


	struct Specie {
		inline static size_t Specie_N = 0;

		Genome repr;
		size_t size = 0;
		size_t num = 0;

		float best_fitness = FLT_MIN;
		size_t gen_since_improv = 0;

		std::vector<size_t> idx_in_population;
	};
	std::vector<Specie> species;


	struct Result {
		Genome* g = nullptr;
		float fitness = 0;
		float adjusted_fitness = 0;
	};
	std::vector<Result> results;

	void complete_with(Genome g, size_t to) noexcept;

	void add_genome(Genome g) noexcept;

	void evaluate(std::function<float(Network&)> fitness) noexcept;
	void evaluate(std::function<float(Network&, void*)> statefull_fitness, void* user) noexcept;
	void speciate() noexcept;
	void select() noexcept;
	void populate() noexcept;

	Genome mutation(const Genome& in) noexcept;

	void turnament_evaluate(

		std::function<float(Network&)> fitness, size_t game_size
	) noexcept;
};
