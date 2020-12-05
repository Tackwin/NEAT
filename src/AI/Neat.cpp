#include "Neat.hpp"
#include <algorithm>
#include <unordered_map>
#include <optional>

double randomd() noexcept { return (double)rand() / (1.0 + RAND_MAX); }
double randomd2() noexcept { auto d = randomd(); return d * d; }

Genome Genome::mate(const Genome& a, const Genome& b) noexcept {
	Genome g;

	size_t i = 0;
	size_t j = 0;
	for (; i < a.connect_genes.size() && j < b.connect_genes.size(); ) {
		auto a_id = a.connect_genes[i].innov;
		auto b_id = b.connect_genes[j].innov;

		if (a_id == b_id) {
			auto r = randomd() < 0.5;
			g.connect_genes.push_back(r ? a.connect_genes[i] : b.connect_genes[j]);
			i++;
			j++;
		} else {
			// disjoint genes are inheret by the fitter parent.
			if (a_id < b_id) {
				g.connect_genes.push_back(a.connect_genes[i]);
				i++;
			}
			if (a_id > b_id) j++;
		}
	}
	for (; i < a.connect_genes.size(); ++i) g.connect_genes.push_back(a.connect_genes[i]);
	for (; j < b.connect_genes.size(); ++j) g.connect_genes.push_back(b.connect_genes[j]);

	std::uint32_t max_nodes = (std::uint32_t)std::max(a.node_genes.size(), b.node_genes.size());

	g.node_genes.resize(max_nodes);

	for (i = 0; i < a.node_genes.size(); ++i) g.node_genes[i] = a.node_genes[i];
	for (; i < b.node_genes.size(); ++i) g.node_genes[i] = b.node_genes[i];

	return g;
}

float relu(float x) noexcept {
	return x > 0 ? x : 0;
}

float lin(float x) noexcept {
	return x;
}

float sig(float x) noexcept {
	return 1 / (1 + std::expf(-x));
}

float activate(Node_Activation_Func f, float x) noexcept {
	switch (f) {
	case Node_Activation_Func::Relu:
		return relu(x);
	case Node_Activation_Func::Lin:
		return lin(x);
	case Node_Activation_Func::Sig:
		return sig(x);
	default:
		return 0;
	}
}

void Network::compute(float* in, size_t in_size, float* out, size_t out_size) noexcept {
	for (size_t i = 0, j = 0; i < nodes.size(); ++i) {
		nodes[i].y = 0;
		if (node_kinds[i] == Node_Kind::Sensor && j < in_size) {
			nodes[i].y = in[j];
			++j;
		}
	}

	for (size_t n = 0; n < 10; ++n) {
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
	for (auto& x : connect_genes) if (x.enabled) if (i-- == 0) x.w = w;
}

Network Genome::phenotype() noexcept {
	Network net;

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
	connect_genes.push_back(c);
}

void Genome::del_connection(size_t i) noexcept {
	connect_genes.erase(std::begin(connect_genes) + i);
}

void Genome::dis_connection(size_t i) noexcept {
	connect_genes[i].enabled = false;
}

float Genome::dist(const Genome& a, const Genome& b, float c1, float c2, float c3) noexcept {
	size_t n_dis = 0;
	size_t n_exc = 0;

	size_t i = 0;
	size_t j = 0;

	float w = 0;
	size_t n_w = 0;
	size_t N = std::max(a.node_genes.size(), b.node_genes.size());

	if (a.connect_genes.empty() && b.connect_genes.empty()) return 0;
	if (a.connect_genes.empty() || b.connect_genes.empty()) return FLT_MAX;

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

		if (i >= a.connect_genes.size()) {
			n_exc = b.connect_genes.size() - j;
			break;
		}
		if (j >= b.connect_genes.size()) {
			n_exc = a.connect_genes.size() - i;
			break;
		}
	}

	if (n_w > 0) w /= n_w;

	return (c1 * n_exc + c2 * n_dis) / N + w;
}

void Neat::complete_with(Genome g, size_t to) noexcept {
	bool need_to_speciate = to != population.size();
	for (size_t i = population.size(); i < to; ++i) add_genome(g);
	if (need_to_speciate) {
		speciate();
	}
}


void Neat::evaluate(std::function<float(Network&)> fitness) noexcept {
	results.clear();

	auto adjust_fitness_of = [&](float x, size_t i) -> float {
		auto n_in_specie = species[genome_info[i].specie].size;

		auto f = x;
		x /= 1 + specie_crowd_rate * n_in_specie;
		x /= 1 + complexity_cost   * population[i].node_genes.size();
		x /= 1 + young_advantage   * genome_info[i].age;
		return x;
	};

	for (size_t i = 0; i < population.size(); ++i) {
		auto& x = population[i];
		auto net = x.phenotype();
		
		Result r;
		r.g = &x;
		r.fitness = fitness(net);
		r.adjusted_fitness = adjust_fitness_of(r.fitness, i);

		results.push_back(r);
	}
}

void Neat::select() noexcept {
	//thread_local std::vector<std::pair<size_t, Result>> indices;
	thread_local std::vector<Genome>      sorted_population;
	thread_local std::vector<Genome_Info> sorted_info;

	//indices.clear();
	sorted_population.clear();
	sorted_info.clear();

	for (auto& x : species) {
		std::sort(
			std::begin(x.idx_in_population),
			std::end(x.idx_in_population),
			[&](size_t a, size_t b) {
				return results[a].adjusted_fitness > results[b].adjusted_fitness;
			}
		);

		x.idx_in_population.resize(x.idx_in_population.size() * survival_rate);
		x.size = x.idx_in_population.size();
	}

	for (auto& x : species) {
		size_t beg = sorted_population.size();

		for (auto& i : x.idx_in_population) {
			sorted_population.push_back(population[i]);
			sorted_info.push_back(genome_info[i]);
		}

		x.idx_in_population.clear();
		for (size_t i = 0; i < x.size; ++i, ++beg) x.idx_in_population.push_back(beg);
	}

	population.swap(sorted_population);
	genome_info.swap(sorted_info);

}


Genome Neat::mutation(const Genome& in) noexcept {
	Genome mutated = in;

	// update connection weight
	if (randomd() < mutation_rate && !mutated.connect_genes.empty()) {
		size_t i = rand() % mutated.connect_genes.size();

		float f = 0.05f * randomd() + 0.975f;
		float b = 0.1f * randomd();

		mutated.update_weight(i, f * mutated.connect_genes[i].w + b);
	}

	// disable connection
	if (5 * randomd() < mutation_rate && !mutated.connect_genes.empty()) {
		size_t i = rand() % mutated.connect_genes.size();
		mutated.dis_connection(i);
	}

	// add connection
	if (5 * randomd() < mutation_rate && mutated.node_genes.size() > 1) {
		std::uint32_t in = rand() % mutated.node_genes.size();
		std::uint32_t out = rand() % (mutated.node_genes.size() - 1);
		float w = (float)(2 * randomd() - 1);
		if (out >= in) out++;

		bool found = false;
		for (auto& x : mutated.connect_genes) if (x.in_id == in && x.out_id == out) {

			found = true;
			break;
		}
		if (!found) if (
			mutated.node_genes[in].kind != Genome::Node_Gene::Output &&
			mutated.node_genes[out].kind != Genome::Node_Gene::Sensor
		) {
			mutated.add_connection(in, out, w);
		}
	}

	// update node function
	if (5 * randomd() < mutation_rate && mutated.node_genes.size() > 0) {
		size_t i = rand() % mutated.node_genes.size();
		Node_Activation_Func f = (Node_Activation_Func)(rand() % Node_Activation_Func::Count);
		mutated.node_genes[i].f = f;
	}

	// add node
	if (15 * randomd() < mutation_rate && !mutated.connect_genes.empty()) {
		// pick a connection
		size_t c = rand() % mutated.connect_genes.size();

		mutated.add_node(c);
	}

	// del node
	if (10 * randomd() < mutation_rate && !mutated.connect_genes.empty()) {
		// pick a connection
		size_t c = rand() % mutated.node_genes.size();

		if (mutated.node_genes[c].kind == Genome::Node_Gene::Kind::Hidden) {
			mutated.del_node(c);
		}

	}


	return mutated;
}

void Neat::populate() noexcept {
	size_t survived = population.size();

	size_t to_give_birth = population_size - population.size();

	for (size_t s_id = 0; s_id < species.size(); ++s_id) if (species[s_id].size > 1) {
		auto& s = species[s_id];

		size_t specie_to_give_birth = (size_t)std::ceil(to_give_birth / (float)species.size());
		specie_to_give_birth = std::min(specie_to_give_birth, to_give_birth);
		to_give_birth -= specie_to_give_birth;

		for (size_t i = 0; i < specie_to_give_birth; ++i) {
			// select parent1 and parent2 from same specie.
			// we are taking a random number squared so the distribution is skewed towards
			// the fittest.
			auto p1_id = (size_t)(randomd2() * s.size);
			auto p2_id = (size_t)(randomd2() * (s.size - 1)); if (p2_id >= p1_id) p2_id++;

			// look for their indices in the population array
			size_t a_idx = species[s_id].idx_in_population[p1_id];
			size_t b_idx = species[s_id].idx_in_population[p2_id];

			// a must be the fittest
			if (results[a_idx].adjusted_fitness < results[b_idx].adjusted_fitness)
				std::swap(a_idx, b_idx);

			auto offspring = Genome::mate(population[a_idx], population[b_idx]);
			offspring = mutation(offspring);

			add_genome(offspring);
		}
	}

	for (size_t i = 0; i < population.size(); ++i) genome_info[i].age++;
}

void Neat::speciate() noexcept {
	for (auto& x : species) x.size = 0;

	std::vector<int> found_specie;
	found_specie.resize(population.size(), -1);

	// First we allocate species based on past representant
	for (size_t s_idx = 0; s_idx < species.size(); ++s_idx) {
		auto& specie = species[s_idx];

		for (size_t i = 0; i < population.size(); ++i) if (found_specie[i] < 0) {
			auto d = Genome::dist(specie.repr, population[i], c_1, c_2, c_3);
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
}

void Neat::add_genome(Genome g) noexcept {
	population.emplace_back(std::move(g));
	Genome_Info i;
	i.age = 0;
	i.specie = 0;
	genome_info.push_back(i);
}

float xor_fitness(Network& net) noexcept {
	float inputs[] = {
		1, 0, 0,
		1, 1, 0,
		1, 0, 1,
		1, 1, 1
	};
	float target_outputs[] = {
		0, 1, 1, 0
	};
	float predicted_outputs[4];

	net.compute(inputs + 0, 3, predicted_outputs + 0, 1);
	net.compute(inputs + 3, 3, predicted_outputs + 1, 1);
	net.compute(inputs + 6, 3, predicted_outputs + 2, 1);
	net.compute(inputs + 9, 3, predicted_outputs + 3, 1);

	float s = 0;
	for (size_t i = 0; i < 4; ++i) {
		float dt = (target_outputs[i] - predicted_outputs[i]);
		s += dt * dt;
	}

	return 1 / (0.001 + s);
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

	g.input_nodes  = *reinterpret_cast<const size_t*>(d + 0);
	g.output_nodes = *reinterpret_cast<const size_t*>(d + 4);

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
