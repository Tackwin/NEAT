#include "Neat.hpp"
#include <algorithm>

double randomd() noexcept {
	return (double)rand() / (double)RAND_MAX;
}

Genome Genome::mate(const Genome& a, const Genome& b) noexcept {
	return {};
}

float relu(float x) noexcept {
	return x > 0 ? x : 0;
}

float activate(Network::Node::Func f, float x) noexcept {
	switch (f) {
	case Network::Node::Func::Relu:
		return relu(x);
	default:
		return 0;
	}
}

void Network::compute(float* in, size_t in_size, float* out, size_t out_size) noexcept {
	for (size_t i = 0, j = 0; i < nodes.size() && j < in_size; ++i) {
		if (node_kinds[i] == Node_Kind::Sensor) {
			nodes[i].y = in[j];
			++j;
		}
	}

	for (size_t n = 0; n < 10; ++n) {
		for (auto& c : connections) nodes[c.out].y += nodes[c.in].y * c.w;

		for (auto& n : nodes) n.y = activate(n.f, n.y);
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
		n.f = n.Relu;

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

void Neat::fill_with(Genome g) noexcept {
	n_inputs = 0;
	n_outputs = 0;
	for (auto& x : g.node_genes) if (x.kind == x.Sensor) n_inputs++;
	for (auto& x : g.node_genes) if (x.kind == x.Output) n_outputs++;

	for (auto& x : population) x = g;
}

// >TODO(Tackwin):
void Neat::add_agents(size_t n) noexcept {}

void Neat::evaluate(std::function<float(Network&)> fitness) noexcept {
	results.clear();

	for (auto& x : population) {
		auto net = x.phenotype();
		
		Result r;
		r.g = &x;
		r.fitness = fitness(net);

		results.push_back(r);
	}
}

void Neat::select() noexcept {
	thread_local std::vector<std::pair<size_t, float>> indices;
	thread_local std::vector<Genome> sorted_population;

	indices.clear();
	sorted_population.clear();

	for (size_t i = 0; i < population.size(); ++i) indices.push_back({i, results[i].fitness});

	std::sort(std::begin(indices), std::end(indices), [] (auto a, auto b) {
		return a.second > b.second;
	});

	indices.resize(population_size * survival_rate);
	for (size_t i = 0; i < indices.size(); ++i)
		sorted_population.push_back(population[indices[i].first]);
	population.swap(sorted_population);
}

void Neat::populate() noexcept {
	size_t survived = population.size();
	for (size_t i = survived; i < population_size; ++i) {
		auto parent = population[i - survived];

		if (randomd() < mutation_rate && !parent.connect_genes.empty()) {
			size_t i = rand() % parent.connect_genes.size();

			float f = 0.05f * randomd() + 0.975f;
			float b = 0.1f * randomd();

			parent.update_weight(i, f * parent.connect_genes[i].w + b);
		}

		if (randomd() < mutation_rate / 5 && !parent.connect_genes.empty()) {
			size_t i = rand() % parent.connect_genes.size();
			parent.dis_connection(i);
		}

		if (randomd() < mutation_rate / 5 && parent.node_genes.size() > 1) {
			std::uint32_t in = rand() % parent.node_genes.size();
			std::uint32_t out = rand() % (parent.node_genes.size() - 1);
			float w = (float)(2 * randomd() - 1);
			if (out >= in) out++;

			if (
				parent.node_genes[in].kind != Genome::Node_Gene::Output &&
				parent.node_genes[out].kind != Genome::Node_Gene::Sensor
			) parent.add_connection(in, out, w);
		}

		if (10 * randomd() < mutation_rate && !parent.connect_genes.empty()) {
			// pick a connection
			size_t c = rand() % parent.connect_genes.size();

			parent.add_node(c);
		}

		population.push_back(parent);
	}
}