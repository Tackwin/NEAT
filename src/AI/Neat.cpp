#include "Neat.hpp"
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include "xstd.hpp"

Genome Genome::mate(const Genome& a, const Genome& b) noexcept {
	Genome g;

	g.serial_number = Genome_Serial++;

	size_t i = 0;
	size_t j = 0;
	for (; i < a.connect_genes.size() && j < b.connect_genes.size(); ) {
		auto a_id = a.connect_genes[i].innov;
		auto b_id = b.connect_genes[j].innov;

		if (a_id == b_id) {
			auto r = randomd() < 0.5;
			g.add_connection(r ? a.connect_genes[i] : b.connect_genes[j]);
			i++;
			j++;
		} else {
			// disjoint genes are inheret by the fitter parent.
			if (a_id < b_id) {
				g.add_connection(a.connect_genes[i]);
				i++;
			}
			if (a_id > b_id) j++;
		}
	}
	for (; i < a.connect_genes.size(); ++i) g.add_connection(a.connect_genes[i]);
	for (; j < b.connect_genes.size(); ++j) g.add_connection(b.connect_genes[j]);

	for (size_t i = 0; i < a.input_nodes; ++i) g.add_new_input();
	for (size_t i = 0; i < a.output_nodes; ++i) g.add_new_output();

	// We construct he minimal set of nodes needed to have every connections.
	std::unordered_set<size_t> hidden_nodes;
	for (auto& x : g.connect_genes) {
		if (x.in_id >= a.input_nodes + a.output_nodes) hidden_nodes.insert(x.in_id);
		if (x.out_id >= a.input_nodes + a.output_nodes) hidden_nodes.insert(x.out_id);
	}

	std::unordered_map<size_t, size_t> nodes_new_idx;
	for (size_t i = 0; i < a.input_nodes + a.output_nodes; ++i) nodes_new_idx[i] = i;

	for (auto& n : hidden_nodes) {
		nodes_new_idx[n] = g.node_genes.size();
		g.add_new_hidden();
	}

	for (auto& x : g.connect_genes) {
		x.in_id = nodes_new_idx[x.in_id];
		x.out_id = nodes_new_idx[x.out_id];
	}

	return g;
}

static float sig(float x)  noexcept { return -1 + 2 / (1 + std::expf(-4.8f * x)); }
static float per(float x)  noexcept { return x > 0.5 ? 1 : 0; }
static float relu(float x) noexcept { return x > 0 ? x : 0; }
static float lin(float x)  noexcept { return x; }

float activate(Node_Activation_Func f, float x) noexcept {
	switch (f) {
	case Node_Activation_Func::Relu:
		return relu(x);
	case Node_Activation_Func::Lin:
		return lin(x);
	case Node_Activation_Func::Sig:
		return sig(x);
	case Node_Activation_Func::Per:
		return per(x);
	default:
		return 0;
	}
}

void Network::clear() noexcept {
	for (auto& x : nodes) x.y = 0;
};

void Network::compute_clear(float* in, size_t in_size, float* out, size_t out_size) noexcept {
	compute(in, in_size, out, out_size);
	clear();
}

void Network::compute(float* in, size_t in_size, float* out, size_t out_size) noexcept {
	for (size_t i = 0, j = 0; i < nodes.size() && j < in_size; ++i) {
		if (node_kinds[i] == Node_Kind::Sensor) {
			nodes[i].y = in[j];
			++j;
		}
	}

	for (size_t n = 0; n < 1; ++n) {
		for (auto& c : connections) nodes[c.out].s += nodes[c.in].y * c.w;
		for (size_t i = 0; i < nodes.size(); ++i) if (node_kinds[i] != Node_Kind::Sensor)
			nodes[i].y = activate(nodes[i].f, nodes[i].s);
		for (auto& n : nodes) n.s = 0;
	}

	for (size_t i = 0, j = 0; i < nodes.size() && j < out_size; ++i) {
		if (node_kinds[i] == Node_Kind::Output) {
			out[j] = nodes[i].y;
			++j;
		}
	}
}

void Genome::add_new_input() noexcept {
	Node_Gene n;
	input_nodes++;
	n.kind = Node_Gene::Kind::Sensor;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

void Genome::add_new_output() noexcept {
	Node_Gene n;
	n.kind = Node_Gene::Kind::Output;
	output_nodes++;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

void Genome::add_new_hidden() noexcept {
	Node_Gene n;
	n.kind = Node_Gene::Kind::Hidden;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

void Genome::add_node(size_t c) noexcept {
	if (!connect_genes[c].enabled) return;
	connect_genes[c].enabled = false;

	Node_Gene n;
	n.kind = n.Hidden;
	n.node_id = node_genes.size();

	node_genes.push_back(n);

	add_connection(connect_genes[c].in_id, n.node_id, connect_genes[c].w / 2);
	add_connection(n.node_id, connect_genes[c].out_id, connect_genes[c].w / 2);
}

void Genome::del_node(size_t c) noexcept {
	node_genes.erase(std::begin(node_genes) + c);
	for (size_t i = connect_genes.size() - 1; i + 1 > 0; --i) {
		if (connect_genes[i].in_id == c || connect_genes[i].out_id == c) {
			connect_genes.erase(std::begin(connect_genes) + i);
		} else {
			if (connect_genes[i].in_id > c) connect_genes[i].in_id--;
			if (connect_genes[i].out_id > c) connect_genes[i].out_id--;
		}
	}
}

void Genome::update_weight(size_t i, float w) noexcept {
	connect_genes[i].w = w;
}

Network Genome::phenotype() const noexcept {
	Network net;
	net.genome_serial_number = serial_number;

	net.nodes.reserve(node_genes.size());
	for (auto& x : node_genes) {
		Network::Node_Kind k;
		Network::Node n;
		n.y = 0;
		n.f = x.f;

		if (x.kind == x.Hidden) k = Network::Node_Kind::Hidden;
		if (x.kind == x.Output) k = Network::Node_Kind::Output;
		if (x.kind == x.Sensor) k = Network::Node_Kind::Sensor;

		net.nodes.push_back(n);
		net.node_kinds.push_back(k);
	}

	net.connections.reserve(connect_genes.size());
	for (auto& x : connect_genes) if (x.enabled) {
		Network::Connections c;
		c.in = x.in_id;
		c.out = x.out_id;
		c.w = x.w;
		net.connections.push_back(c);
	}

	return net;
}

void Genome::add_connection(std::uint32_t in, std::uint32_t out, float w) noexcept {
	Connect_Gene c;
	c.in_id = in;
	c.out_id = out;
	c.w = w;
	c.innov = Innov_N++;

	add_connection(c);
}

void Genome::add_connection(const Connect_Gene& c) noexcept {
	auto it = std::find_if(
		std::begin(connect_genes),
		std::end(connect_genes),
		[&](const auto& a) { return a.in_id == c.in_id && a.out_id == c.out_id; }
	);

	if (it != std::end(connect_genes)) it->w += c.w;
	else connect_genes.push_back(c);
}


void Genome::del_connection(size_t i) noexcept {
	connect_genes.erase(std::begin(connect_genes) + i);
}

void Genome::dis_connection(size_t i) noexcept {
	connect_genes[i].enabled = false;
}

float Genome::dist(const Genome& a, const Genome& b, float c1, float c2, float c3) noexcept {
	if (a.connect_genes.empty() && b.connect_genes.empty()) return 0;
	if (a.connect_genes.empty() || b.connect_genes.empty()) return FLT_MAX;

	if (a.connect_genes.back().innov < b.connect_genes.back().innov)
		return dist(b, a, c1, c2, c3);

	size_t n_dis = 0;
	size_t n_exc = 0;

	size_t i = 0;
	size_t j = 0;

	float w = 0;
	size_t n_w = 0;
	size_t N = std::max(a.node_genes.size(), b.node_genes.size());
	if (N <= 20) N = 1;
	else        N = N - 19;


	for (; i < a.connect_genes.size() && j < b.connect_genes.size(); ) {
		auto a_id = a.connect_genes[i].innov;
		auto b_id = b.connect_genes[j].innov;

		if (a_id == b_id) {
			auto d = a.connect_genes[i].w - b.connect_genes[j].w;
			w += d > 0 ? d : -d;

			++i;
			++j;
			n_w++;
		}
		if (a_id > b_id) {
			n_dis++;
			j++;
		}
		if (b_id > a_id) {
			n_dis++;
			i++;
		}
	}

	n_exc = a.connect_genes.size() - std::min(i, a.connect_genes.size());

	if (n_w > 0) w /= n_w;

	return (c1 * n_exc + c2 * n_dis) / N + c3 * w;
}

void Neat::complete_with(Genome g, size_t to) noexcept {
	bool need_to_speciate = to != population.size();
	for (size_t i = population.size(); i < to; ++i) add_genome(g);
	if (need_to_speciate) {
		speciate();
	}
}


void Neat::evaluate(std::function<float(Network&)> fitness) noexcept {
	evaluate([fitness](Network& a, void*) { return fitness(a); }, nullptr);
}

void Neat::evaluate(std::function<float(Network&, void*)> statefull_fitness, void* user) noexcept {
	results.clear();

	for (size_t i = 0; i < population.size(); ++i) {
		auto& x = population[i];
		auto net = x.phenotype();
		
		Result r;
		r.g = &x;
		r.fitness = statefull_fitness(net, user);

		results.push_back(r);
	}
}

void Neat::select() noexcept {
	thread_local std::vector<size_t>      indices;
	thread_local std::vector<Genome_Info> sorted_genome_info;
	thread_local std::vector<Genome>      sorted_population;
	thread_local std::vector<Result>      sorted_results;

	sorted_genome_info.clear();
	sorted_population.clear();
	sorted_results.clear();
	indices.clear();

	size_t total_kill = population_size * (1 - survival_rate);
	size_t population_wide_kill = total_kill * population_competition_rate;
	size_t specie_wide_kill = total_kill - population_wide_kill;

	// species kill
	for (auto& specie : species) {

		std::sort(
			std::begin(specie.idx_in_population),
			std::end(specie.idx_in_population),
			[&] (size_t a, size_t b) { return results[a].fitness > results[b].fitness; }
		);

		size_t this_specie_kill = std::ceil(specie.size * (1 - population_competition_rate));
		this_specie_kill = std::min(this_specie_kill, specie_wide_kill);
		this_specie_kill = std::min(this_specie_kill, specie.size - 1);

		specie.idx_in_population.resize(specie.size - this_specie_kill);

		for (auto& i : specie.idx_in_population) {
			sorted_genome_info.emplace_back(std::move(genome_info[i]));
			sorted_population.emplace_back(std::move(population[i]));
			sorted_results.emplace_back(std::move(results[i]));

			i = sorted_population.size() - 1;
		}
	}

	for (size_t i = 0; i < sorted_population.size(); ++i) indices.push_back(i);
	std::sort(
		std::begin(indices),
		std::end(indices),
		[&] (size_t a, size_t b) {
			return sorted_results[a].adjusted_fitness > sorted_results[b].adjusted_fitness;
		}
	);

	size_t cutoff_idx = indices[indices.size() - population_wide_kill - 1];
	indices.resize(indices.size() - population_wide_kill);

	for (auto& specie : species) {
		
		for (size_t i = specie.idx_in_population.size() - 1; i + 1 > 0; --i)
			if (specie.idx_in_population[i] > cutoff_idx)
				specie.idx_in_population.erase(std::begin(specie.idx_in_population) + i);

		specie.size = specie.idx_in_population.size();
		if (specie.size > 0) specie.repr = population[specie.idx_in_population.front()];
	}

	genome_info.clear();
	population.clear();
	results.clear();

	for (auto& i : indices) {
		genome_info.emplace_back(std::move(sorted_genome_info[i]));
		population.emplace_back(std::move(sorted_population[i]));
		results.emplace_back(std::move(sorted_results[i]));
	}



	/*
	indices.clear();

	for (auto& x : species) {
		std::sort(
			std::begin(x.idx_in_population),
			std::end(x.idx_in_population),
			[&](size_t a, size_t b) {
				return results[a].adjusted_fitness > results[b].adjusted_fitness;
			}
		);

		size_t to_kill = x.idx_in_population.size() * (1 - survival_rate);
		to_kill *= (1 - population_competition_rate);

		x.idx_in_population.resize(x.idx_in_population.size() - to_kill);
		x.size = x.idx_in_population.size();
	}

	for (auto& x : species) {
		size_t beg = sorted_population.size();

		for (auto& i : x.idx_in_population) {
			sorted_population.emplace_back(std::move(population[i]));
			sorted_info.emplace_back(std::move(genome_info[i]));
			sorted_results.push_back(results[i]);
		}

		x.idx_in_population.clear();
		for (size_t i = 0; i < x.size; ++i, ++beg) x.idx_in_population.push_back(beg);
	}

	for (size_t i = 0; i < sorted_population.size(); ++i) indices.push_back(i);

	std::sort(
		std::begin(indices),
		std::end(indices),
		[&](size_t a, size_t b) {
			return sorted_results[a].adjusted_fitness > sorted_results[b].adjusted_fitness;
		}
	);

	population.clear();
	genome_info.clear();

	size_t to_kill = sorted_population.size() * (1 - survival_rate);
	to_kill *= population_competition_rate;

	for (size_t i = 0; i + to_kill < sorted_population.size(); ++i) {
		population.emplace_back(std::move(sorted_population[indices[i]]));
		genome_info.emplace_back(std::move(sorted_info[indices[i]]));
	}

	for (auto& x : species) x.idx_in_population.clear();
	for (size_t i = 0; i < population.size(); ++i)
		species[genome_info[i].specie].idx_in_population.push_back(i);
	for (auto& x : species) x.size = x.idx_in_population.size();

	*/
}


Genome Neat::mutation(const Genome& in) noexcept {
	Genome mutated = in;

	auto& p = params;

	// update connection weight
	for (auto& x : mutated.connect_genes) {
		if (randomd() < p.update_connection_rate * mutation_rate) {
			if (randomd() < p.new_connection_weight_rate * mutation_rate)
				x.w = randomd() * 2 - 1;
			else
				x.w += randomd2() * 2 - 1;
		}
	}

	for (auto& x : mutated.connect_genes)
		if (randomd() < p.toggle_connection_rate * mutation_rate) x.enabled = !x.enabled;

	// add connection
	if (randomd() < p.add_connection_rate * mutation_rate && mutated.node_genes.size() > 1) {
		std::uint32_t in  = (std::uint32_t)(randomd() * mutated.node_genes.size());
		std::uint32_t out = (std::uint32_t)(randomd() * mutated.node_genes.size());
		float w = (float)(2 * randomd2() - 1);
		
		// no if (out >= in) out++; here because we might want recurrent connections.

		if (
			mutated.node_genes[in].kind != Genome::Node_Gene::Output &&
			mutated.node_genes[out].kind != Genome::Node_Gene::Sensor
		) {
			mutated.add_connection(in, out, w);
		}
	}

	for (auto& x: mutated.node_genes) if (randomd() < p.update_node_act_rate * mutation_rate)
		x.f = (Node_Activation_Func)(rand() % Node_Activation_Func::Count);

	// add node
	if (randomd() < p.add_node_rate * mutation_rate && !mutated.connect_genes.empty()) {
		// pick a connection
		size_t c = rand() % mutated.connect_genes.size();

		mutated.add_node(c);
	}

	for (size_t i = mutated.node_genes.size() - 1; i + 1 > 0; --i) {
		if (randomd() < p.del_node_rate * mutation_rate) {
			mutated.del_node(i);
		}
	}

	return mutated;
}

void Neat::populate() noexcept {
	size_t survived = population.size();

	size_t to_give_birth = population_size - population.size();
#if 0
	for (size_t s_id = 0; s_id < species.size(); ++s_id) {
		auto& s = species[s_id];

		size_t specie_to_give_birth = (size_t)std::ceil(to_give_birth / (float)species.size());
		specie_to_give_birth = std::min(specie_to_give_birth, to_give_birth);
		to_give_birth -= specie_to_give_birth;

		if (s.size > 1) {
			for (size_t i = 0; i < specie_to_give_birth; ++i) {
				// select parent1 and parent2 from same specie.
				// we are taking a random number squared so the distribution is skewed towards
				// the fittest.
				auto p1_id = (size_t)(randomd2() * s.size);
				auto p2_id = (size_t)(randomd2() * (s.size - 1)); if (p2_id >= p1_id) p2_id++;

				// look for their indices in the population array
				size_t a_idx = s.idx_in_population[p1_id];
				size_t b_idx = s.idx_in_population[p2_id];

				// a must be the fittest
				if (results[a_idx].adjusted_fitness < results[b_idx].adjusted_fitness)
					std::swap(a_idx, b_idx);

				auto offspring = Genome::mate(population[a_idx], population[b_idx]);
				offspring = mutation(offspring);

				add_genome(offspring);
			}
		} else {
			for (size_t i = 0; i < specie_to_give_birth; ++i) {
				add_genome(mutation(population[s.idx_in_population.front()]));
			}
		}
	}
#endif

	for (size_t i = 0; i < to_give_birth; ++i) {
		auto p1_id = (size_t)(randomd2() * population.size());
		auto p2_id = (size_t)(randomd2() * (population.size() - 1)); if (p2_id >= p1_id) p2_id++;

		auto* a = &population[p1_id];
		auto* b = &population[p2_id];
		
		if (results[p1_id].adjusted_fitness < results[p2_id].adjusted_fitness) std::swap(a, b);

		add_genome(mutation(Genome::mate(*a, *b)));
	}

	for (size_t i = 0; i < population.size(); ++i) genome_info[i].age++;
	for (size_t i = 0; i < population.size(); ++i) results[i].g = &population[i];
}

void Neat::speciate() noexcept {
	for (auto& x : species) x.size = 0;

	std::vector<int> found_specie;
	found_specie.resize(population.size(), -1);

	// First we allocate species based on past representant
	for (size_t s_idx = 0; s_idx < species.size(); ++s_idx) {
		auto& specie = species[s_idx];

		for (size_t i = 0; i < population.size(); ++i) if (found_specie[i] < 0) {
			float d = Genome::dist(specie.repr, population[i], c_1, c_2, c_3);
			if (d > specie_dt) continue;

			found_specie[i] = s_idx;
			species[s_idx].size++;
			genome_info[i].specie = s_idx;
		}
	}

	size_t new_specie_idx = species.size();

	// After that if there is still people who have not been assigned a specie.
	for (size_t i = 0; i < population.size(); ++i) if (found_specie[i] < 0) {
		size_t s_idx = 0;

		size_t j = new_specie_idx;
		for (; j < species.size(); ++j) {
			auto& specie = species[j];

			auto d = Genome::dist(specie.repr, population[i], c_1, c_2, c_3);
			if (d > specie_dt) continue;

			s_idx = j;
			break;
		}

		if (j == species.size()) {
			Specie new_specie;
			new_specie.size = 0;
			new_specie.repr = population[i];
			new_specie.num = Specie::Specie_N++;

			s_idx = species.size();
			species.push_back(new_specie);
		}

		found_specie[i] = s_idx;
		species[s_idx].size++;
		genome_info[i].specie = s_idx;
	}

	for (size_t i = 0; i < species.size(); ++i) if (species[i].size == 0) {
		species.erase(std::begin(species) + i);
		for (auto& x : genome_info) if (x.specie >= i) x.specie--;

		--i;
	}

	for (auto& x : species) x.idx_in_population.clear();
	for (size_t i = 0; i < population.size(); ++i)
		species[genome_info[i].specie].idx_in_population.push_back(i);

	if (specifie_number_of_species) {
		// Here we can adjust specie_dt if the user has set a preferred number of species.

		if (preferred_number_of_species > species.size())      specie_dt /= 1.5;
		else if (preferred_number_of_species < species.size()) specie_dt *= 1.5;
	}
	
	// adjust fitness
	auto adjust_fitness_of = [&](float x, size_t i) -> float {
		auto n_in_specie = species[genome_info[i].specie].size;
		n_in_specie = std::log2(n_in_specie + min_specie_size_advantage);

		size_t complexity = population[i].node_genes.size();

		auto f = x;
		x /= 1 + complexity_cost   * complexity;
		x /= 1 + specie_crowd_rate * n_in_specie;
		x /= 1 + young_advantage   * genome_info[i].age;
		return x;
	};
	for (size_t i = 0; i < population.size(); ++i)
		results[i].adjusted_fitness = adjust_fitness_of(results[i].fitness, i);
}

void Neat::add_genome(Genome g) noexcept {
	population.emplace_back(std::move(g));
	Genome_Info i;
	i.age = 0;
	i.specie = 0;
	genome_info.push_back(i);

	Result r;
	r.adjusted_fitness = 0;
	r.fitness = 0;
	r.g = &population.back();
	results.push_back(r);
}

std::vector<std::uint8_t> Genome::serialize() const noexcept {
	std::vector<std::uint8_t> data;

	std::uint32_t input  = input_nodes;
	std::uint32_t output = output_nodes;

	data.resize(16);
	memcpy(data.data() + 0, &input, 4);
	memcpy(data.data() + 4, &output, 4);

	std::uint32_t n_nodes = node_genes.size();
	std::uint32_t n_connect = connect_genes.size();

	memcpy(data.data() + 8, &n_nodes, 4);
	memcpy(data.data() + 12, &n_connect, 4);

	for (auto& x : node_genes) {
		size_t i = data.size();
		data.resize(i + sizeof(Genome::Node_Gene));
		memcpy(data.data() + i, &x, sizeof(Genome::Node_Gene));
	}

	for (auto& x : connect_genes) {
		size_t i = data.size();
		data.resize(i + sizeof(Genome::Connect_Gene));
		memcpy(data.data() + i, &x, sizeof(Genome::Connect_Gene));
	}

	return data;
}

Genome Genome::deserialize(const std::vector<std::uint8_t>& data) noexcept {
	auto* d = data.data();
	Genome g;

	g.input_nodes  = *reinterpret_cast<const std::uint32_t*>(d + 0);
	g.output_nodes = *reinterpret_cast<const std::uint32_t*>(d + 4);

	g.node_genes.resize(*reinterpret_cast<const std::uint32_t*>(d + 8));
	g.connect_genes.resize(*reinterpret_cast<const std::uint32_t*>(d + 12));


	for (size_t i = 0; i < g.node_genes.size(); ++i) {
		g.node_genes[i] =
			*reinterpret_cast<const Genome::Node_Gene*>(d + 16 + i * sizeof(Genome::Node_Gene));
	}
	
	for (size_t i = 0; i < g.connect_genes.size(); ++i) {
		g.connect_genes[i] = *reinterpret_cast<const Genome::Connect_Gene*>(
			d + 16
				+ g.node_genes.size() * sizeof(Genome::Node_Gene)
				+ i * sizeof(Genome::Connect_Gene)
		);
	}

	return g;
}
