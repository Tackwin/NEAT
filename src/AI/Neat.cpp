#include "Neat.hpp"


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
	for (size_t i = 0, j = 0; j < in_size && i < input_nodes.size(); ++i) {
		nodes[input_nodes[i]].y = in[j];
	}

	for (size_t n = 0; n < 10; ++n) {
		for (auto& c : connections) nodes[c.out].y += nodes[c.in].y * c.w;

		for (auto& n : nodes) n.y = activate(n.f, n.y);
	}

	for (size_t i = 0, j = 0; i < output_nodes.size() && i < out_size; ++i) {
		out[j] = nodes[output_nodes[i]].y;
	}
}

void Genome::add_new_input() noexcept {
	Node_Gene n;
	n.kind = Node_Gene::Kind::Sensor;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

void Genome::add_new_output() noexcept {
	Node_Gene n;
	n.kind = Node_Gene::Kind::Output;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

void Genome::add_node() noexcept {
	Node_Gene n;
	n.kind = n.Hidden;
	n.node_id = node_genes.size();
	node_genes.push_back(n);
}

Network Genome::phenotype() noexcept {
	Network net;

	net.nodes.reserve(node_genes.size());
	for (auto& x : node_genes) {
		Network::Node n;
		n.y = 0;
		n.f = n.Relu;

		if (x.kind == x.Sensor) net.input_nodes.push_back(net.nodes.size());
		if (x.kind == x.Output) net.output_nodes.push_back(net.nodes.size());

		net.nodes.push_back(n);
	}

	net.connections.reserve(connect_genes.size());
	for (auto& x : connect_genes) if (x.enabled) {
		Network::Connections c;
		c.in = x.in_id;
		c.out = x.out_id;
		c.w = c.w;
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
